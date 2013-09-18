#! /bin/bash
#|
exec mzscheme -qu "$0" ${1+"$@"}
|#
;;====================================================================================================
;; This is a script to parse .dnd tree files.

;; NOTE: This was the "Version 1" script and has been replaced by the haskell version, phybin.hs

;;====================================================================================================

;; [2009.12.06] Hacking it to work with wolbachia_A_group/input

(module find_tree_matches mzscheme  
  ;; Import the parser and lexer generators.
  (require (lib "yacc.ss" "parser-tools")
	   (lib "lex.ss" "parser-tools")
	   (prefix : (lib "lex-sre.ss" "parser-tools"))
	   (lib "list.ss")
	   ;(lib "iu-match.ss")
	   "iu-match.ss"
	   )
;  (provide (all-defined))

(require (lib "pretty.ss")
	 (lib "process.ss")
	 (prefix srfi: (lib "1.ss" "srfi")))
  
;;====================================================================================================
;; First define the PARSER
;;====================================================================================================

(define-empty-tokens op-tokens 
   (LeftParen RightParen Comma SemiColon Colon
    EOF ))

(define-tokens value-tokens (NAME NUM))

(define-lex-abbrevs (lower-letter (:/ #\a #\z))
                    (upper-letter (:/ #\A #\Z))
		    (alpha (:or lower-letter upper-letter))
		    (digit (:/ #\0 #\9))
                    (digits (:+ digit))
		    (letter (:or alpha digit "." "'" "-" "/"))
		    
		    ;; This requires that names have an underscore in them:
		    ;; That's what distinguishes them from numbers:
		    ;(name (:seq (:+ letter) (:+ (:seq "_" (:* letter)))))
		    
		    ;; [2010.08.13] I can't remember why I did that... were there gene names that started with numbers??
		    (name (:seq alpha (:+ (:or letter "_"))))
		    
		    ;; [2009.12.06] Odd, getting negative branch lengths!!
		    ;; [2009.12.07] Ran into some trees with 0 distance and no decimal point:
		    (number (:seq (:? "-") digits (:? (:seq "." digits (:? (:seq "E-" digits))))))
		    ;(number (:seq digit digits))
		    )

(define lex
  (lexer-src-pos

   [(eof) 'EOF]
   
   ;; Ignore all whitespace:
   [(:or #\tab #\space #\newline #\return) (return-without-pos (lex input-port))]

   ["(" 'LeftParen]
   [")" 'RightParen]  
   ["," 'Comma] 
   [";" 'SemiColon] 
   [":" 'Colon]

   ;; Like most of these scripts this is a hack.  It assumes all node ID's in the tree start with a particular string.
   
   ;[(:seq "c" name) (token-NAME lexeme)]
   [name (token-NAME lexeme)]

   [number (token-NUM (string->number lexeme))]
   ; [number (token-NUM  lexeme)]
   
   ))

(define (format-pos pos)
  (if (position-line pos)
      (format "line ~a:~a" (position-line pos) (position-col pos))
      (format "char ~a" (position-offset pos))))

(define parse
 (parser        
   (src-pos)
   (start wholefile)
   (end EOF)

  (tokens value-tokens op-tokens)
  (error (lambda (a b c start end) 
            (fprintf (current-error-port)
		    "PARSE ERROR: after ~a token, ~a~a.\n" 
                    (if a "valid" "invalid") b 
                    (if c (format " carrying value ~s" c) ""))
            (fprintf (current-error-port)
		     "  Located between ~a and ~a.\n"
                    (format-pos start) (format-pos end))))
   ;; Precedence:
  ;(precs  	  )
   ;(debug "_parser.log")

;    Tree: The full input Newick Format for a single tree
;    Subtree: an internal node (and its descendants) or a leaf node
;    Leaf: a leaf node
;    Internal: an internal node (and its descendants)
;    BranchSet: a set of one or more Branches
;    Branch: a tree edge and its descendant subtree.
;    Name: the name of a node
;    Length: the length of a tree edge.
; [edit]The grammar rules

(grammar

   (wholefile [(Subtree SemiColon) $1]
	      [(Branch SemiColon)  $1])
   (Subtree [(Leaf) $1]
            [(Internal) $1])
   (Leaf [(Name) $1])
   (Internal [(LeftParen BranchSet RightParen Name) 
	      (if $4 (cons $4 $2) $2)])

   (BranchSet [(Branch) (list $1)]
              ;[(BranchSet Comma Branch) (cons $1 $3)]
	      [(Branch Comma BranchSet) (cons $1 $3)]
	      )
   (Branch [(Subtree Length) (list $2 $1)])

   (Name [() #f]
         [(NAME) $1])
   (Length [() 0]
           [(Colon NUM) $2])
)

#;
   (grammar

    (wholefile [(tree SemiColon) $1]
	       ;; ACK making semicolon and final tag optional!!
	       [(tree) $1]
	       [(LeftParen nodes+ RightParen) (cons 0.0 $2)]
	       [(LeftParen nodes+ RightParen SemiColon) (cons 0.0 $2)]
	       )

    (numtag [(Colon NUM) $2]
            ;; [2010.08.13] Making them optional everywhere: 
            [() 0.0]
	    )


    (tree [(NAME numtag) (list $2 $1)]
	  [(LeftParen nodes+ RightParen numtag) (cons $4 $2)])
    
    (nodes+ [(tree) (list $1)]
	    [(tree Comma nodes+) (cons $1 $3)])
)

))

;; ================================================================================

(define (parse-file f) 
    (let ([p (open-input-file f)])    
      (let ([res (parse-port p)])
        ;(printf "PARSED: \n")(pretty-print  res)
	(close-input-port p)
	res)))

(define (parse-port ip)
  (port-count-lines! ip)
  (parse (lambda () 
            (let ((lexed (lex ip)))
	      ;(printf "LEXED: ~s\n" (position-token-token lexed))
	      (flatten lexed))))
  )

(define (flatten pt)
  (let loop ((x pt))
        (if (position-token? (position-token-token x))
            (begin (error 'flatten "No nested position-tokens!")
              (loop (position-token-token x)))
            x)))

(define allargs (vector->list (current-command-line-arguments)))
(define pretty? #t)
 ;; End main

(print-graph #t)
(print-vector-length #f)


;;====================================================================================================
;; Now normalize the trees.
;;====================================================================================================

;; The basic idea here is to sort branches by some criteria.
;; If that criteria turns out to be insufficient... then we can add extra criteria.
(define (normalize-tree tr)
  ;; This returns a branch-sorted tree
  (define (loop tr )
    (match tr
      [(,dist ,name) (guard (string? name))
        (values (list dist name) 1 1)]
      [(,dist ,x)
       (error 'normalize-tree "was not expecting interior nodes with one sub-branch... ~s" (list dist x))]
      [(,dist ,[nodes sizes depths] ...)
       ;; This creates an ugly data type (list overuse) with both the branches and their size/depth measurement:
       (define sorted 
	 (sort 
	  (map list nodes sizes depths)
	  (lambda (ls1 ls2) 
	    ;(printf "Comparing ~s and ~s\n" (cdr ls1) (cdr ls2))
	    (or (< (cadr ls1) (cadr ls2))
		(and (= (cadr ls1) (cadr ls2))
		     (< (caddr ls1) (caddr ls2)))))
	  ))	
        ;(printf "\n ~a Node sizes ~s and depths ~s\n" dist sizes depths)
	;; If there are any non-leaf "plateaus" in the ramp then we can't deal with that yet:
       (let loop ((last (car sorted)) 
		  (rest (cdr sorted)))
	 (cond 
	  [(null? rest) (void)]
	  ;; If they're equal according to the size/depth metric:
	  [(and (equal? (cdr last) (cdar rest))
		;; And non-leaf:
		(not (equal? (cdr last) '(1 1))))
	   ;; Ok, given that we have a match, that's still ok if they are the SAME tree. 
	   ;; The only problem is if they're different in a way that our metric doesn't capture.
	   (if (tree-shape-equals? (car last) (caar rest))
	       (void)
	       (begin
	         (printf "\nERRORful subbranches:\n")
		 (pretty-print (car last))
		 (pretty-print (caar rest))
		 (error 'normalize-tree "Cannot yet handle same depth/size signatures of sub-branches: ~s ~s" 
			(cdr last) (cdar rest))))]
	  [else (loop (car rest) (cdr rest))]))

        (values (cons dist (map car sorted))
		(apply + 1 sizes)
		(add1 (apply max depths))
		)]))
  (define-values (tree _ __) (loop tr))
  tree)

;; Assumes that the trees are already in normal form:
;; Compares tree shape and leaf nodes, IGNORES interior nodes.
(define (tree-equals? tr1 tr2)
  ;(printf "EQUALS? ~s ~s\n" tr1 tr2)
  (match (cons tr1 tr2)
    [((,dist1 ,name1) . (,dist2 ,name2))
     (guard (string? name1) (string? name2))
     (if 
          #t                    ;; -- This would give you only shape.
	  (string=? name1 name2) ;; -- This requires that the leaf names be the same.
	  )
     ]
    [((,dist1 ,nodes1 ...) . (,dist2 ,nodes2 ...))
     (and (equal? (length nodes1) (length nodes2))
	  (andmap tree-equals? nodes1 nodes2))]
    [(,s1 . ,s2) (guard (string? s1) (string? s2)) (error 'tree-equals? "implementation assumption violated... oops")]
    [,else #f]))

(define (tree-shape-equals? tr1 tr2)
  ;(printf "\nSHAPE EQUALS? ~s ~s\n" tr1 tr2)
  (match (cons tr1 tr2)
    [((,dist1 ,name1) . (,dist2 ,name2))
     (guard (string? name1) (string? name2))
     #t ;; -- This would give you only shape.
     ]
    [((,dist1 ,nodes1 ...) . (,dist2 ,nodes2 ...))
     (and (equal? (length nodes1) (length nodes2))
	  (andmap tree-shape-equals? nodes1 nodes2))]
    [(,s1 . ,s2) (guard (string? s1) (string? s2)) (error 'tree-shape-equals? "implementation assumption violated... oops")]
    [,else #f]))


;;====================================================================================================
;; TESTS
;;====================================================================================================

(define tests 
 '(
   "(,,(,));"                               ; no nodes are named
   "(A,B,(C,D));"                           ; leaf nodes are named
   "(A,B,(C,D)E)F;"                         ; all nodes are named
   "(:0.1,:0.2,(:0.3,:0.4):0.5);"           ; all but root node have a distance to parent
   "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;"       ; all have a distance to parent
   "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"       ; distances and leaf names (popular)
   "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"     ; distances and all names
   "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;"    ; a tree rooted on a leaf node (rare)
   ))

(define (test)  (map parse tests))

;;====================================================================================================
;; Here's our script "main"
;;====================================================================================================

(define (new-main files)
  (define __ (begin (printf "\nBeginning processing of ~s files\n" (length files))))
  ; (define parsed (parse-file file))
  ; (define normalized (normalize-tree parsed))
  ;(pretty-print parsed) ;(newline)
  ;(printf "NORMALIZED ~s\n " file)
  ;(pretty-print normalized)
  ;;(exit)

  (define ___ (begin (printf "Parsing and normalizing trees: ")(flush-output)))
  (define all-trees
    (map (lambda (file)
           ;(define __ (begin (printf " ~a " file)(flush-output)))
	   (define parsed (parse-file file))
	   (define normalized (normalize-tree parsed))
	   (printf ".") (flush-output)
	   (list file normalized))
      files))

  ;; Type: (((name tree) ...) ...)
  ;; Each equivalence class is represented by (tree name ...)
  ;; Only a single representative tree is needed.
  (define classes '())

  (if (null? (cdr all-trees))
      (begin
        ;; If we have one input just print it:
	(printf "\nOnly one file as input.  Not making equivalence classes, just printing it:\n")
	(pretty-print (cadar all-trees))        
	)
      (begin 
	(printf "\nComputing equivalence classes on ~s trees: " (length all-trees)) (flush-output)
	(for-each
	    (lambda (entry)
	      (define found-hit #f)
	      (match entry
		[(,file ,normtree)
		 (let ([new (map (lambda (class) 
				   (if (tree-equals? normtree (car class))
				       (begin 
					 (set! found-hit #t)
					 (cons (car class) (cons file (cdr class))))
				       class))
			      classes)])
		   (printf ".") (flush-output)		    
		   (set! classes
			 (if found-hit new
			     (cons `(,normtree ,file) classes))))
		 ]))
	  all-trees)
	(newline)

	;; Now sort them:
	(set! classes (sort classes (lambda (a b) (< (length a) (length b)))))

	(printf "Finished.  ~s files fell into ~s equivalence classes\n" (length files) (length classes))
	(printf "Sizes of each equivalence class:  ")
	(for-each (lambda (class) (printf "  ~s" (length (cdr class)))) classes)
	(newline)
	(printf "Dumping lists of files to class<N>.txt, and representative tree to class<N>.tree...\n")

	(system "rm -f class*.txt class*.tree")
	
	(for-each (lambda (i class) 
		;(when (file-exists? (format "class~s.txt" i))  (delete-file (format "class~s.txt" i)))
		;(when (file-exists? (format "class~s.tree" i)) (delete-file (format "class~s.tree" i)))
		(with-output-to-file (format "class~s.txt" i)
		  (lambda ()
		    (for-each (lambda (file) (printf "~a " file)) (cdr class))
		    (newline)))
		(with-output-to-file (format "class~s.tree" i)
		  (lambda () (pretty-print (car class))))
		)
      (srfi:iota (length classes))
      classes)
)))

(define (print-usage)
  (printf "Usage: find_tree_matches <file-to-process> ...\n\n"))


;COUNTING: ((0 ((0 "lpg140") (0 "lpl135"))) (0 "lpp135") (0 "LPC_081"))

(define (sumcount stuff) 
  ;(printf "sumcount ~s\n" stuff)
  (if (or (string? stuff) (not stuff)) 0
      (foldl (lambda (x acc) (+ acc (count x))) 0 stuff)))

(define (count tree)
  ;(define __ (printf "COUNTING: ~s\n" tree))
     (if (number? (car tree))
	  (+ 1 (sumcount (cadr tree)))
	  (sumcount tree)))

;; BENCHMARK TEST:

;(for-each (lambda (file) (parse-file file)) allargs)

(pretty-print (foldl + 0 (map (lambda (file) (count (parse-file file))) allargs)))
;(pretty-print (map (lambda (file) (pretty-print (count (parse-file file)))) allargs))

;(pretty-print (count '((0 ((0 "lpg292") (0 "lpl285"))) (0 "lpp299") (0 "LPC_323"))))



#;
(begin
  (when (null? allargs) 
    (print-usage) (flush-output)
  ;(error 'find_tree_matches "no files given")
  (printf "Error: no files given!\n")
  (exit 1))
  (new-main allargs))


) ;; End module.


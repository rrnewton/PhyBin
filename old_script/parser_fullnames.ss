#! /bin/bash
#|
exec mzscheme -qu "$0" ${1+"$@"}
|#

;; I'm duplicating the code from parser.ss and hacking it.

(module parser_fullnames mzscheme  
  ;; Import the parser and lexer generators.
  (require (lib "yacc.ss" "parser-tools")
	   (lib "lex.ss" "parser-tools")
	   (prefix : (lib "lex-sre.ss" "parser-tools"))
	   (lib "list.ss")
	   ;(lib "iu-match.ss")
	   "iu-match.ss"
	   )
;  (provide (all-defined))

(require (lib "pretty.ss"))
  
(define-empty-tokens op-tokens 
   (LeftParen RightParen Comma SemiColon Colon
    EOF ))

(define-tokens value-tokens (NAME NUM))

(define-lex-abbrevs (lower-letter (:/ "a" "z"))
                    (upper-letter (:/ #\A #\Z))
		    (digit (:/ "0" "9"))
                    (digits (:* digit))
		    (letter (:or lower-letter upper-letter digit "." "'" "-" "/"))
		    (name (:seq (:+ letter) (:* (:seq "_" (:* letter)))))
		    (number (:seq digits "." digits (:? (:seq "E-" digits))))
		    )

(define ws-lex
  (lexer-src-pos
   [(eof) 'EOF]
   
   ;; Ignore all whitespace:
   [(:or #\tab #\space #\newline #\return) (return-without-pos (ws-lex input-port))]

   ["(" 'LeftParen]
   [")" 'RightParen]  
   ["," 'Comma] [";" 'SemiColon] [":" 'Colon]

   [(:seq "'silva_" name "'") 
    (begin ;(printf "TEST ~a\n" (substring lexeme 7 (sub1 (string-length lexeme))))
	   (token-NAME (substring lexeme 7 (sub1 (string-length lexeme))))
	   )]
   
   [(:seq "silva_" name) (token-NAME (substring lexeme 6 (string-length lexeme)))]


   ;[number (token-NUM (string->number lexeme))]
   [number (token-NUM  lexeme)]
   
   ))

(define (format-pos pos)
  (if (position-line pos)
      (format "line ~a:~a" (position-line pos) (position-col pos))
      (format "char ~a" (position-offset pos))))

(define ws-parse
 (parser        
   (src-pos)
   (start wholefile)
   (end EOF)

  (tokens value-tokens op-tokens)
  (error (lambda (a b c start end) 
            (printf "PARSE ERROR: after ~a token, ~a~a.\n" 
                    (if a "valid" "invalid") b 
                    (if c (format " carrying value ~s" c) ""))
            (printf "  Located between ~a and ~a.\n"
                    (format-pos start) (format-pos end))))
   ;; Precedence:
  ;(precs  	  )
   ;(debug "_parser.log")

   (grammar

    (wholefile [(tree SemiColon) $1]
	       ;; ACK making semicolon and final tag optional!!
	       [(tree) $1]
	       [(LeftParen nodes+ RightParen) (cons "0.0" $2)]
	       )

    (numtag [(Colon NUM) $2])

    (tree [(NAME numtag) (list $2 $1)]
	  [(LeftParen nodes+ RightParen numtag) (cons $4 $2)])
    
    (nodes+ [(tree) (list $1)]
	    [(tree Comma nodes+) (cons $1 $3)])

)))

;; ================================================================================

(define (ws-parse-file f) 
    (let ([p (open-input-file f)])    
      (let ([res (ws-parse-port p)])
	(close-input-port p)
	res)))

(define (ws-parse-port ip)
  (port-count-lines! ip)
  ;; cdr
  (ws-parse (lambda () (flatten (ws-lex ip))))
  ;(position-token-token (ws-lex ip))
  )

(define (flatten pt)
  ;(printf " ")
  (let loop ((x pt))
        (if (position-token? (position-token-token x))
            (begin (error 'flatten "No nested position-tokens!")
              (loop (position-token-token x)))
            x)))

(define allargs (vector->list (current-command-line-arguments)))
(define pretty? #t)
(define (main filename)
  #;
  (if pretty?      
      (pretty-print (ws-parse-file filename))
      (begin (write (ws-parse-file filename))
	     (newline)
	     ))

  (define table '())
  (define counter 1)
  (define (clean-name n)
    (define matches (regexp-match #rx".*__(.*)" n))
    (define cleaned (car (reverse matches)))
    ;(printf "MATCHES ~a\n" matches)
    (set! table `((,counter ,cleaned) . ,table))
    (set! counter (+ 1 counter))
    ;(sub1 counter)

    ;; THIS IS ALL I CHANGED:
    cleaned
    )

  (with-output-to-file (string-append filename ".nameonly")
      (lambda ()	
	(let loop ((x (ws-parse-file filename)))
	  (match x
	    ;; A leaf:
	    [(,dist ,name) (guard (string? name))       
	     (printf "~a:~a" (clean-name name) dist)]
	    [(,dist ,children ...)
	     (printf "(")
	     (let inner ((ls children))
	       (cond 
		[(null? ls)       (void)]
		[(null? (cdr ls)) (loop (car ls))]
		[else (loop (car ls)) (printf ",")
		      (inner (cdr ls))]))
	     (printf "):~a" dist)]
	    ))
	(printf ";\n")))
  
  (with-output-to-file (string-append filename ".table")
    (lambda ()
      (for-each 
	  (lambda (ls) (printf "~a ~a\n" (car ls) (cadr ls)))
	(reverse table))))
  )

(print-graph #t)
(print-vector-length #f)


;; Here's our script invocation:

;; When run in --persist mode we run in a loop.
(if (null? allargs)
	(error 'wsparse "No filename provided...")
	(main (car allargs)))
(printf "Done writing output files.\n")
;(exit 0)

) ;; End module.


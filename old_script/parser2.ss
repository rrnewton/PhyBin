#! /bin/bash
#|
exec chez --script "$0" ${1+"$@"}
|#

(printf "Called with ~s\n" (command-line-arguments))

(include "match.ss")
(import iu-match)

;[2001.07.15]
(define port->slist
  (lambda (p)
    (let porttoslistloop ([exp (read p)] [acc '()])
        (if (eof-object? exp)
            (begin (close-input-port p)
                   (reverse! acc))
            (porttoslistloop (read p) (cons exp acc))))))
;; [2008.05.07] Adding a hack to tolerate #! lines at the top of the file.
(define file->slist
  (lambda (filename . opts)
    (let* ([p (apply open-input-file filename opts)]
	   [line1 (get-line p)])
      ;; Doesn't allow whitespace!
      (if (and (>= (string-length line1) 2)
	       (string=? "#!" (substring line1 0 2)))
	  ;; Read the rest of it straight away:
	  (port->slist p)
	  ;; Oops, we need to "put back" that first line  We do that by just starting over.
	  (begin (close-input-port p)
		 (port->slist (apply open-input-file filename opts))))
      )))

(define sexps (file->slist (car (command-line-arguments))))

(define wstrings
  (match sexps
    [() '()]
    [(,uq ,[x]) (guard (eq? uq 'unquote)) x]
    [,s (guard (symbol? s)) (symbol->string s)]
    [(,[hd] . ,[tl]) `(,hd . ,tl)]
    ))

(pretty-print wstrings)


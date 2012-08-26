;;;; reset_spacetime.lisp
;;;; Authors: Jonah Miller (jonah.maxwell.miller@gmail.com)
;;;;          Rajesh Kommu 
;;;;          David Kamensky

;;;; This file contains functions which reset the simulation. These
;;;; are useful for debugging and for batch operations on spacetimes.

;;;--------------------------------------------------------------------------
;;; RESET THE SIMULATION
;;;--------------------------------------------------------------------------

;;; JM: Resetting the simulation at runtime is tricky business. The
;;; following two functions reset the system in different ways. It's
;;; not terribly important that these functions work properly. They're
;;; mostly usefull for debuging at runtime.

;; reset-spacetime-fast is a faster way to reset the spacetime, but it
;; is dangerous because I am not certain all hash tables and important
;; variables are reset. Use at your own risk.
(defun reset-spacetime-fast ()
  "the state of the simulation after (reset-spacetime) is identical 
to the state after (load \"cdt2p1.lisp\"). Use at your own risk."
  ;; clear the hash tables
  (clrhash *TL2SIMPLEX->ID*)
  (clrhash *SL2SIMPLEX->ID*)
  (clrhash *TL1SIMPLEX->ID*)
  (clrhash *SL1SIMPLEX->ID*)
  (clrhash *ID->SPATIAL-2SIMPLEX*)
  (clrhash *ID->3SIMPLEX*)
  ;; reset the counters
  ;; (setf *LAST-USED-2SXID* 0)
  (setf *LAST-USED-3SXID* 0)
  (setf *RECYCLED-3SX-IDS* '())
  (setf *LAST-USED-POINT* 0)
  ;; reset the ''bulk'' variables
  (setf N0 0)
  (setf N1-SL 0)
  (setf N1-TL 0)
  (setf N2-SL 0)
  (setf N2-TL 0)
  (setf N3-TL-31 0)
  (setf N3-TL-22 0)
  ;; Reset the ''boundary'' variables
  (setf *N1-SL-TOP* 0)
  (setf *N1-SL-BOT* 0)
  (setf *N3-22-TOP* 0)
  (setf *N3-22-BOT* 0)
  (setf *N3-31-TOP* 0)
  (setf *N3-31-BOT* 0)
  ;; reset the parameters
  (setf *k0* 0.0)
  (setf *k3* 0.0)
  (setf *eps* 0.02)
  (setf *a* 1.0)
  (setf *alpha* -1.0)
  (setf *k* 1.0)
  (setf *litL* 1.0)
  (setf CURRENT-MOVE-IDENTIFIER "UNKNOWN")
  (setf CURRENT-MOVE-NUMBER 0)
  (setf STOPOLOGY "unknown")
  (setf BCTYPE "unknown")
  (setf SAVE-EVERY-N-SWEEPS 10)
  (setf NUM-T 666666)
  (setf N-INIT 0)
  (setf NUM-SWEEPS 0)
  ;; reset the move markers
  (setf 26MARKER 0.0)
  (setf 23MARKER 0.0)
  (setf 44MARKER 0.0)
  (setf 32MARKER 0.0)
  (setf 62MARKER 0.0)
  ;; Reset attempted moves list
  (reset-move-counts))

;; reset-spacetime-slow simply reloads all modules and essentially
;; restarts the simulation. Only useful for runtime debugging.
(defun reset-spacetime-slow nil
  "Reloads all modules and restarts simulation. Use at your own risk."
  (load "cdt2p1.lisp"))


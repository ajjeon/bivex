(ns bivex.core
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules]))


(defn -main
  "main function to run bivex"
  []
  ()
  )

; apply changes & move head

(defn change-chrom
  "apply a rule and update the nucleosome, head leaves"
  [rule nuc]
  (assoc nuc
         (keyword (:class rule)) (:right rule)
         :head 0))

(defn drop-nuc
  [no_idx chromtape]
  (remove #(= no_idx %) chromtape))

(defn get-the-rest-idx
  "get the index of the rest of the unchanging nucleosome"
  [prevnuc_idx nextnuc_idx chromtape]
  (drop-nuc nextnuc_idx (drop-nuc prevnuc_idx (range (count chromtape)))))

(defn move-head
  "assign new head position"
  [nuc_idx]
  (assoc (second (nth chromatin/chromtape nextnuc_idx)) :head 1))

; we need to keep both the position of the current head and the next head 
; BE CAREFUL - next head is chosen at random
(defn apply-rule
  [chromtape]
  (let [prevnuc_idx (chromatin/find-idx-with-head chromtape)
      nextnuc_idx (chromatin/nucleo-idx-next-head chromtape)
      prevnuc (chromatin/find-nucleosome-with-head chromtape)
      rule (rules/select-rule prevnuc)
      prevnuc_new [prevnuc_idx (change-chrom rule prevnuc)] 
      nextnuc_new [nextnuc_idx (move-head nextnuc_idx)]
      new_chromtape (concat
                     (map #(nth chromtape %)
                          (get-the-rest-idx prevnuc_idx nextnuc_idx chromtape))
                     (vector nextnuc_new)
                     (vector prevnuc_new))]
    (println prevnuc_idx nextnuc_idx prevnuc rule prevnuc_new nextnuc_new)
  (sort new_chromtape)))



(defn copy-chrom
  "copy the marks from input chromatin"
  [chromatin cidx]
;  (nth (clojure.string/split chromtape #"_") cidx)
  )

;; work on each nucleosome. if applying rule, 


(defn change-chrom
  "for changing nucleosome, change. not, copy over"
  [chrom]
  (let [chrom chrom]
      (cond (= (:head chrom) 1)


            (print "change")
            :else (print "no"))))

(map #(change-chrom (second (nth chromtape %)))
     (take 10 (range)))


;; (if idx == changing index, change the nucleosome
;;    not, copy the nucleosome

(doseq [[idx item] chromtape]
;  (println item)
  (cond (= (:head item) "1") (change-chrom item idx)
  )
)


;; (clojure.walk/prewalk-replace [(second (nth chromtape 2)) [:head "1"]]) ;; updating values

;; (map #(vector % (inc (% chromtape))) [:a :b])

;; (vector 1 (update-in (second (nth chromtape 1)) [:head] inc))

;;;;; at each iteration, recruit-change-move-eval


(cond (= state "RECRUIT") () ;; recruit matching modifiers
      (= state "APPLY") () ;; apply changes
      (= state "MOVE") () ;; move head
      (= state "EVALUATE") () ;; evaluate - iteration, methyl proportion, gene exp outcome
      )



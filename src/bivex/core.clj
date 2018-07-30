(ns bivex.core
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules])
  (:gen-class))


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
  [nextnuc_idx]
  (assoc (second (nth chromatin/chromtape nextnuc_idx)) :head 1))


(defn print-nucleosome-withouthead
  [nucleosome]
  (clojure.string/join "|" [(:k4 (second nucleosome)) (:k27 (second nucleosome))]))

(defn print-nucleosome-withhead
  [nucleosome]
  (clojure.string/join "*" [(:k4 (second nucleosome)) (:k27 (second nucleosome))]))

(defn print-nucleosome
  [nucleosome]
  (cond (= (:head (second nucleosome)) 1) (print-nucleosome-withhead nucleosome)
        :else (print-nucleosome-withouthead nucleosome)))


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
    (println prevnuc_idx nextnuc_idx)
    (println (clojure.string/join "__" (map print-nucleosome (sort new_chromtape))))
    (sort new_chromtape)
    ))


(defn -main
  "main function to run bivex"
  []
  (take 50 (iterate apply-rule chromatin/chromtape))
  )

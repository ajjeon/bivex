(ns bivex.core
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules])
  (:require [jutsu.core :as j])
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
  [chromtape nextnuc_idx]
  (assoc (second (nth chromtape nextnuc_idx)) :head 1))


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

(defn plot-line
  [y]
  (j/graph! "K4 Chart"
  [{:x (range (count y))
    :y y
    :type "scatter"}]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; we need to keep both the position of the current head and the next head 
; BE CAREFUL - next head is chosen at random
(defn apply-rule
  [chromtape]
  (let [prevnuc_idx (chromatin/find-idx-with-head chromtape)
        nextnuc_idx (chromatin/nucleo-idx-next-head chromtape)
        prevnuc (chromatin/find-nucleosome-with-head chromtape)
        rule (rules/select-rule prevnuc)
        prevnuc_new [prevnuc_idx (change-chrom rule prevnuc)] 
        nextnuc_new [nextnuc_idx (move-head chromtape nextnuc_idx)]
        new_chromtape (concat
                       (map #(nth chromtape %)
                            (get-the-rest-idx prevnuc_idx nextnuc_idx chromtape))
                       (vector nextnuc_new)
                       (vector prevnuc_new))]
    (sort new_chromtape)))

(defn evaluate-chrom
  [new_chromtape y]
  (let [y (conj [y] (apply + (map #(:k4 (second %)) new_chromtape)))]
    (plot-line y)
    y
))

(for [iter (range 10)]
  )

(loop [i 0]  
  (when (< i 5)    
;    (let [chromtape (apply-rule chromtape)
    (recur (inc i)); loop i will take this value
))

(defn run
  [chromtape]
  (let [new_chromtape (apply-rule chromtape)
        y (conj [] (apply + (map #(:k4 (second %)) chromtape)))]
;    (println (clojure.string/join "__" (map print-nucleosome (sort new_chromtape))))
    (evaluate-chrom new_chromtape y)
    (println y)
    new_chromtape
    ))


;;;; TODO trying to pass both y value and chromtape
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn -main
  "main function to run bivex"
  []
  (take 50 (iterate apply-rule chromatin/chromtape))
  )

(take 500 (iterate apply-rule chromatin/chromtape))


(defn calculate-ratio

  []
  )


(let [k27 (apply + (map #(:k27 (second %)) chromatin/chromtape))
      k4 (apply + (map #(:k4 (second %)) chromatin/chromtape))]
  (cond (and (not= k27 0) (not= k4 0)) (/ k4 k27)
        (and (= k27 0) (not= k4 0)) 
        (and (= k4 0) (not= k27 0)) (
        (and (= k27 0) (= k4 0)) (0)


        )))

k27 0 k4 0
k27 1 k4 0
k27 0 k4 1
k27 1 k4 1


(plot-line [1 2 3 4 5 6] [8 4 5 8 2 1])



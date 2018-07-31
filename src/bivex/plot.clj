(ns bivex.plot
  (:require [jutsu.core :as j]))

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
  [k4 k27 both genex]
  (j/graph! "Valency Chart"
  [{:x (range (count k4))
    :y k4
    :type "scatter"
    :name "H3K4me3 mono"}
   {:x (range (count k4))
    :y k27
    :type "scatter"
    :name "H3K27me3 mono"}
   {:x (range (count k4))
    :y both
    :type "scatter"
    :name "Bivalent"}
   {:x (range (count k4))
    :y genex
    :type "bar"
    :name "Gene expression"}
   ]))

(defn plot-bar
  [chromtape]
  (j/graph! "Snapshop of chromatin"
            [{:x (range (count chromtape))
              :y (vec (map #(:k4 (second %)) chromtape))
              :type "bar"
              :name "H3K4me3"}
             {:x (range (count chromtape))
              :y (vec (map #(:k27 (second %)) chromtape))
              :type "bar"
              :name "H3K27me3"}]))

(defn plot-cell-genex
  [genex ncells]
  (j/graph! ""
            [{:x [0 1]
              :y genex
              :type "bar"}]))

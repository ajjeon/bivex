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
            [{:x ["OFF" "ON"]
              :y [(second (first (frequencies genex)))
                  (second (second (frequencies genex)))
                  ]
              :type "bar"}]))

(defn plot-cell-all
  [k4mono k27mono biv genex]
  (j/graph! "Snapshop of all cells"
            [{:x (range (count k4mono))
              :y k4mono
              :type "bar"
              :name "H3K4me3 mono"}
             {:x (range (count k4mono))
              :y k27mono
              :type "bar"
              :name "H3K27me3 mono"}
             {:x (range (count k4mono))
              :y biv
              :type "bar"
              :name "Bivalent"}
             {:x (range (count k4mono))
              :y genex
              :type "bar"
              :name "Gene exp"}]))

(defn plot-cell-sum
  [k4mono k27mono biv genex]
  (j/graph! "Bulk reads (sum)"
            [{:x [0]
              :y [(apply + k4mono)]
              :type "bar"
              :name "H3K4me3 mono"}
             {:x [0]
              :y [(apply + k27mono)]
              :type "bar"
              :name "H3K27me3 mono"}
             {:x [0]
              :y [(apply + biv)]
              :type "bar"
              :name "Bivalent"}
             {:x [0]
              :y [(apply + genex)]
              :type "bar"
              :name "Gene exp"}]))

(defn get-sum-marks-nucleosome
[chromtape_array n mark]
(apply + (map #(mark (second (nth (:chromtape %) n))) chromtape_array))
)

(defn plot-bar-sum
  [chromtape_array]
  (let [nnuc (count (:chromtape (nth chromtape_array 1)))
        k4nuc (map #(get-sum-marks-nucleosome chromtape_array % :k4) (range nnuc))
        k27nuc (map #(get-sum-marks-nucleosome chromtape_array % :k27) (range nnuc))]
  (j/graph! "Cumulative histone marks across nucleosomes"
            [{:x (range nnuc)
              :y k4nuc
              :type "scatter"
              :name "H3K4me3"}
             {:x (range nnuc)
              :y k27nuc
              :type "scatter"
              :name "H3K27me3"}])))



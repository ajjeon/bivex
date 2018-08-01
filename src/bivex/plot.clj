(ns bivex.plot
  (:require [jutsu.core :as j]))

(defn print-nucleosome-withouthead
  "prints nucleosomes that are not currently being read"
  [nucleosome]
  (clojure.string/join "|" [(:k4 (second nucleosome)) (:k27 (second nucleosome))]))

(defn print-nucleosome-withhead
  "prints nucleosomes that are currently being read"
  [nucleosome]
  (clojure.string/join "*" [(:k4 (second nucleosome)) (:k27 (second nucleosome))]))

(defn print-nucleosome
  "print the nucleosome based on the head position"
  [nucleosome]
  (cond (= (:head (second nucleosome)) 1) (print-nucleosome-withhead nucleosome)
        :else (print-nucleosome-withouthead nucleosome)))

(defn plot-line
  "plots the traceable valency chart"
  [k4 k27 both genex]
  (j/graph! "Valency Chart"
  [{:x (range (count k4))
    :y k4
    :mode "lines"
    :type "scatter"
    :name "H3K4me3 mono"}
   {:x (range (count k4))
    :y k27
    :mode "lines"
    :type "scatter"
    :name "H3K27me3 mono"}
   {:x (range (count k4))
    :y both
    :mode "lines"
    :type "scatter"
    :name "Bivalent"}
   {:x (range (count k4))
    :y (map #(+ % 9) genex)
    :mode "lines"
    :type "scatter"
    :line {:shape "hvh"}
    :name "Gene expression"}
   ]))

(defn plot-bar
  "graphical representation of the histone mark pattern at one timepoint"
  [chromtape]
  (j/graph! "Snapshot of chromatin"
            [{:x (range (count chromtape))
              :y (vec (map #(:k4 (second %)) chromtape))
              :type "bar"
              :name "H3K4me3"}
             {:x (range (count chromtape))
              :y (vec (map #(:k27 (second %)) chromtape))
              :type "bar"
              :name "H3K27me3"}]))

(defn plot-cell-genex
  "bar plot for the number of cells that are active and silence"
  [genex ncells]
  (j/graph! ""
            [{:x ["OFF" "ON"]
              :y [(second (first (frequencies genex)))
                  (second (second (frequencies genex)))
                  ]
              :type "bar"}]))

(defn plot-cell-all
  "bar plot to show both histone mark patterns and the gene expression status across cells."
  [k4mono k27mono biv genex]
  (j/graph! "Snapshot of all cells"
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
              :y (map #(+ % 9) genex)
              :mode "lines+markers"
              :type "scatter"
              :line {:shape "hvh"}
              :name "Gene exp"}]))

(defn plot-cell-sum
  "histone marks are normalized by the total number of nucleosome (i.e. ncells * nnucleosomes, default nnucleosomes is 10. genex value is normalized by the total number of cells"
  [k4mono k27mono biv genex ncells]
  (j/graph! "Bulk reads (sum)"
            [{:x [0]
              :y [(format "%3f" (* (float (/ (apply + k4mono) (* ncells 10))) 100))]
              :type "bar"
              :name "H3K4me3 mono"}
             {:x [0]
              :y [(format "%3f" (* (float (/ (apply + k27mono) (* ncells 10))) 100))]
              :type "bar"
              :name "H3K27me3 mono"}
             {:x [0]
              :y [(format "%3f" (* (float (/ (apply + biv) (* ncells 10))) 100))]
              :type "bar"
              :name "Bivalent"}
             {:x [0]
              :y [(format "%3f" (* (float (/ (apply + genex) ncells)) 100))]
              :type "bar"
              :name "Gene exp"}]))

(defn get-sum-marks-nucleosome
  "calculates how many 'mark' is present at each chromtape"
  [chrom_in n mark]
  (apply + (map #(mark (second (nth (:chromtape %) n))) chrom_in))
)

(defn plot-bar-sum
  "scatter plot for cumulative histone marks across nucleosomes from many cells"
  [chrom_in]
  (let [nnuc (count (:chromtape (nth chrom_in 1)))
        k4nuc (map #(get-sum-marks-nucleosome chrom_in % :k4) (range nnuc))
        k27nuc (map #(get-sum-marks-nucleosome chrom_in % :k27) (range nnuc))]
    (j/graph! "Cumulative histone marks across nucleosomes"
              [{:x (range nnuc)
                :y k4nuc
                :mode "lines+markers"
                :type "scatter"
                :name "H3K4me3"}
               {:x (range nnuc)
                :y k27nuc
                :mode "lines+markers"
                :type "scatter"
                :name "H3K27me3"}])))



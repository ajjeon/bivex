(ns bivex.plot
  (:require [jutsu.core :as j])
  (:require [clojure.string :as string]))

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

(defn print-chromtape
  [chromtape]
  (println (clojure.string/join "__" (map print-nucleosome (sort chromtape)))))

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
    :y (map #(cond (zero? %) 11 :else 12) genex)
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
              :y (map #(cond (zero? %) 11 :else 12) genex)
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
              :y [(format "%3f" (* (float (/ (apply + (map #(cond (zero? %) 1 :else 0) genex)) ncells)) 100))]
              :type "bar"
              :name "Gene exp"}]))

(defn get-sum-marks-nucleosome
  "calculates how many 'mark' is present at each chromtape"
  [chrom_in n mark]
  (apply + (map #(mark (second (nth (:chromtape %) n))) chrom_in)))

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

;; when-switch function should first see if before "afteriter" the genex value was zero. if yes, check when the switch to on happened after "afteriter". if not, zero
;;; if switch to on never happened after "afteriter", zero.

(defn when-switch?*
  "find the last occurence of ON gene expression and add 1 to indicate when the stable switch to activation happened"
  [cvector beforeiter] 
  (cond (reduce identical? cvector) (+ beforeiter 1) 
        :else (+ (string/index-of (string/join (map str cvector)) "11111") (+ beforeiter 1))) ; find the first occurence of "stable" -- at least 5 consecutive iters -- reactivation

  )

(defn when-switch?
  [genex beforeiter]
  (let [sub-genex (vec (map #(cond (zero? %) 0 :else 1) genex))
        sub-genex-pre (subvec sub-genex beforeiter)]
    ;; (println sub-genex)
    ;; (println sub-genex-pre)
    (cond (and (= (nth sub-genex beforeiter) 0) (= (last sub-genex) 1))
          (cond (some #(= % 1) sub-genex-pre) (when-switch?* sub-genex-pre beforeiter)
          :else nil)
          :else nil)))

(defn when-switch-plot
  [ncells allgenex beforeiter]
  (let [switches (vec (map #(when-switch? (:genex %) beforeiter) allgenex))]
    (println switches)
    (j/graph!
     "When stable switch happened in each cell"
     [{:x (range ncells)
       :y switches
       :mode "markers"
       :type "scatter"
       :name "switch"}
      {:x (range ncells)
       :y (repeat ncells beforeiter) 
       :mode "lines"
       :type "scatter"
       :name "EZH2i added"}])))


(defn plot-nucleo-mat
  [save-chromtape geneidx chromtape_count genex]
;  (println genex)
  (j/graph!
   "Mark evolution"
   [{:x [0]
     :z save-chromtape
     :type "heatmap"
     ;:colorscale "Jet"
     }
    {:x [(first geneidx) (inc (last geneidx))]
     :y [-100 -100]
     :mode "lines"
     :type "scatter"}
    {:x [-2]
     :z (apply mapv vector [(vec (map #(cond (zero? %) 0 :else 1) genex))])
     :type "heatmap"
     :showscale false
     :colorscale "Viridis"}]
   {:autosize "false"
     :width 1200
     :height 1200}))


(ns bivex.eval
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules])
  (:require [bivex.plot :as plot])
  (:require [bivex.cell :as cell])
  (:require [jutsu.core :as j])
  (:gen-class))


(defn evaluate-chrom-one
  "evaluate chromtape for one cell"
  [chrom_in]
  (let [new_chrom_in (cell/apply-rule chrom_in)
        new_chromtape (:chromtape new_chrom_in)
        y (into (:k4mono chrom_in) (:k4mono (cell/check-valency new_chromtape)))
        y2 (into (:k27mono chrom_in) (:k27mono (cell/check-valency new_chromtape)))
        y3 (into (:biv chrom_in) (:biv (cell/check-valency new_chromtape)))
        genex (into (:genex chrom_in) (cell/gene-on? (last y) (last y2) (last y3)))]
    (plot/plot-line y y2 y3 genex) ; trace valency and gene expression
    (plot/plot-bar new_chromtape) ; snapshot of valency
    (Thread/sleep 100)
    (println (clojure.string/join "__" (map plot/print-nucleosome (sort new_chromtape)))) ; trace iteration
;    (println new_chrom_in)
    {:k4mono y
     :k27mono y2
     :biv y3
     :genex genex
     :chromtape new_chromtape
     :rules (:rules new_chrom_in)}  ;;; UPDATE RULES HERE
))

(defn evaluate-chrom-bulk
  "evaluate chromtape multiple cells"
  [chrom_in]
  (let [new_chrom_in (cell/apply-rule chrom_in)
        new_chromtape (:chromtape new_chrom_in)
        y (into (:k4mono chrom_in) (:k4mono (cell/check-valency new_chromtape)))
        y2 (into (:k27mono chrom_in) (:k27mono (cell/check-valency new_chromtape)))
        y3 (into (:biv chrom_in) (:biv (cell/check-valency new_chromtape)))
        genex (into (:genex chrom_in) (cell/gene-on? (last y) (last y2) (last y3)))]
    {:k4mono y
     :k27mono y2
     :biv y3
     :genex genex
     :chromtape new_chromtape
     :rules (:rules new_chrom_in)}  ;;; UPDATE RULES HERE
))

;;;; cells iterate in parallel. capture snapshot after each iteration
(defn run-one
  "run one cell through iterations"
  [chrom_in itern]
  (last (take itern (iterate evaluate-chrom-one chrom_in)))
)

(defn run-bulk
  "run multiple cells through iterations"
  [chrom_in ncells itern]
  (let [cell_group (repeat ncells chrom_in)
        t (pmap #(last (take itern (iterate evaluate-chrom-bulk %))) cell_group)
        x (map #(select-keys % [:k4mono :k27mono :biv :genex]) t)
        t2   (map #(map last (vals %)) x )
        k4mono (vec (map #(nth % 0) (vec t2)))
        k27mono (vec (map #(nth % 1) (vec t2)))
        biv (vec (map #(nth % 2) (vec t2)))
        genex (vec (map #(nth % 3) (vec t2)))]
;    (println (apply + k4mono) (apply + k27mono) (apply + biv) (apply + genex))
    (plot/plot-cell-all k4mono k27mono biv genex)
    (plot/plot-cell-sum k4mono k27mono biv genex ncells)
    (plot/plot-bar-sum t)
))


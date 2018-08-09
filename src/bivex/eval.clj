(ns bivex.eval
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules])
  (:require [bivex.plot :as plot])
  (:require [bivex.cell :as cell])
  (:require [bivex.files :as files])
  (:require [jutsu.core :as j])
  (:gen-class))

(defn generate_chrom_in [rfile cfile trate]
  (let [chromtape (files/read-in-chromatin cfile)
        k4mono (:k4mono (cell/check-valency chromtape))
        k27mono (:k27mono (cell/check-valency chromtape))
        biv (:biv (cell/check-valency chromtape))
        rules (files/read-in-file rfile)
        ]
;    (println trate)
    {:k4mono k4mono
     :k27mono k27mono
     :biv biv
     :genex (cell/gene-on? (first k4mono) (first k27mono) (first biv) trate)
     :chromtape chromtape
     :trate trate
     :rules rules
     :orules rules}
    ))

(defn evaluate-chrom
  "applies a selected rule and update the chrom_in"
  [chrom_in]
  (let [new_chrom_in (cell/apply-rule chrom_in)
        new_chromtape (:chromtape new_chrom_in)
        y (into (:k4mono chrom_in) (:k4mono (cell/check-valency new_chromtape)))
        y2 (into (:k27mono chrom_in) (:k27mono (cell/check-valency new_chromtape)))
        y3 (into (:biv chrom_in) (:biv (cell/check-valency new_chromtape)))
        genex (into (:genex chrom_in) (cell/gene-on? (last y) (last y2) (last y3) (:trate chrom_in)))]
    {:k4mono y
     :k27mono y2
     :biv y3
     :genex genex
     :chromtape new_chromtape
     :rules (:rules new_chrom_in)
     :orules (:orules chrom_in)
     :trate (:trate chrom_in)}
))

(defn evaluate-chrom-with-plot
  "plots updated chrom_in. only applies to single-cell iterations"
  [chrom_in]
  (let [new_chrom_in (evaluate-chrom chrom_in)]
  (plot/plot-line (:k4mono new_chrom_in) (:k27mono new_chrom_in) (:biv new_chrom_in) (:genex new_chrom_in)) ; trace valency and gene expression
;   (plot/plot-bar new_chromtape) ; snapshot of valency
  (Thread/sleep 100)
;  (println (clojure.string/join "__" (map plot/print-nucleosome (sort (:chromtape new_chrom_in))))) ; trace iteration
                                        ;    (println new_new_chrom_in)
  new_chrom_in
    ))

;;;; for single-cell simulations
(defn run-one
  "run one cell through iterations"
  [chrom_in itern]
  (last (take itern (iterate evaluate-chrom-with-plot chrom_in)))
)

(defn run-one-with-change
  [init_chrom_in beforeiter afteriter]
  (let [new-rules (files/read-in-file @rules/new-rules-file)
        beforerun (run-one init_chrom_in beforeiter)]
;      (reset! rules/default-rules-file @rules/new-rules-file)
    (run-one (assoc
              (assoc beforerun :rules new-rules)
              :orules new-rules) afteriter)
      ))

;;;; for multiple-cell simulations

(defn run-many
  "run one cell through iterations"
  [chrom_in itern]
  (last (take itern (iterate evaluate-chrom chrom_in)))
)

(defn run-many-with-change
  [init_chrom_in beforeiter afteriter]
  (let [new-rules (files/read-in-file @rules/new-rules-file)
        beforerun (run-many init_chrom_in beforeiter)]
;      (reset! rules/default-rules-file @rules/new-rules-file)
    (run-many (assoc
               (assoc beforerun :rules new-rules)
               :orules new-rules) afteriter)
      ))

(defn run-bulk
  "run multiple cells through iterations"
  [chrom_in ncells beforeiter afteriter]
  (let [cell_group (map #(generate_chrom_in @rules/default-rules-file @chromatin/chromatin-file %) (vec (repeatedly ncells #(+ (rand-int 9) 1))))
        t (pmap #(run-many-with-change % beforeiter afteriter) cell_group)
        x (map #(select-keys % [:k4mono :k27mono :biv :genex]) t)
        t2 (map #(map last (vals %)) x)
        allgenex (map #(select-keys % [:genex]) t)
        k4mono (vec (map #(nth % 0) (vec t2)))
        k27mono (vec (map #(nth % 1) (vec t2)))
        biv (vec (map #(nth % 2) (vec t2)))
        genex (vec (map #(nth % 3) (vec t2)))
        ]
;    (println allgenex)
;    (println (apply + k4mono) (apply + k27mono) (apply + biv) (apply + genex))
    ;; (plot/plot-cell-all k4mono k27mono biv genex)
    ;; (plot/plot-cell-sum k4mono k27mono biv genex ncells)
    ;; (plot/plot-bar-sum t)
    ;; (plot/when-switch-plot ncells allgenex beforeiter)
    (println (map #(last (:genex %)) allgenex)) 
;    (plot/genex-box-plot allgenex)
    ))


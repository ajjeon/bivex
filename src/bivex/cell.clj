(ns bivex.cell
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules])
  (:require [bivex.tools :as tools]))

 ;(cond (= (:action rule) "methyltransferase") 1 :else 0)

;; (defn move-head
;;   [chromtape nuc_h]
;;   (let [nuc_n_idx (chromatin/nucleo-idx-next-head chromtape (first nuc_h))
;;         nuc_n_new [nuc_n_idx (tools/cast-head chromtape nuc_n_idx)]
;;         nuc_rest (tools/chromtape-rest chromtape [nuc_n_idx])]
;;     (sort (concat [nuc_n_new] nuc_rest))))

(defn apply-rule
  [chromtape rule nuc_h]
  (let [nuc_h_new [(first nuc_h) (tools/change-chrom rule (second nuc_h))]
        nuc_rest (tools/chromtape-rest chromtape [(first nuc_h)])]
    (sort (concat [nuc_h_new] nuc_rest))))

;; (defn apply-rule
;;   [chrom_in]
;;   ;; find nucleosome with head, select a rule
;;   (let [nuc_h (chromatin/find-nucleosome-with-head (:chromtape chrom_in))
;;         rule (rules/select-rule (second nuc_h) (:rules chrom_in))
;;         eval_chromtape (cond (= (:action rule) "turnover")
;;                              (turnover (:chromtape chrom_in) rule nuc_h)
;;                              :else
;;                              (apply-rule-sub (:chromtape chrom_in) rule nuc_h))
;;         new_chromtape (move-head eval_chromtape nuc_h)
;; ;        new_rule (rules/update-rules new_chromtape (:orules chrom_in) nuc_h)
;;         ]
;; ;    (println rule)
;;   {:k4mono (:k4mono chrom_in)
;;    :k27mono (:k27mono chrom_in)
;;    :biv (:biv chrom_in)
;;    :genex (:genex chrom_in)
;;    :chromtape (sort new_chromtape)
;;    :rules (vec new_rule) 
;;    :orules (:orules chrom_in)}))

(defn check-valency
  "check valency of the given chromtape"
  [new_chromtape]
  (let [k4total (vec (map #(:k4 (second %)) new_chromtape))
        k27total (vec (map #(:k27 (second %)) new_chromtape))
        biv (map #(bit-and %1 %2) k4total k27total)
        xor (map #(bit-xor %1 %2) k4total k27total)
        k4mono (map #(bit-and %1 %2) k4total xor)
        k27mono (map #(bit-and %1 %2) k27total xor)]
    {:k4mono (vector (apply + k4mono))
     :k27mono (vector (apply + k27mono))
     :biv (vector (apply + biv))}))

(defn gene-on?
     "calls gene expression outcome and multiplies by the transcriptional rate"
     [k4mono k27mono biv trate]
  (cond (> k4mono (+ k27mono 1)) [trate] ;; this is where we set the gene expression threshold
        :else [0]))

(defn call-valency-per-cell
  "call valency state per cell, based on the largest value."
  [chrom_in]
  (let [k4 {:k4mono (last (:k4mono chrom_in))}
        k27 {:k27mono (last (:k27mono chrom_in))} 
        biv {:biv (last (:biv chrom_in))}            
        genex (last (:genex chrom_in))
        markmap (merge k4 k27 biv)
        maxmark (cond (every? zero? (vals markmap)) "NA"
                      (not (apply distinct? (vals markmap))) "NA"
                      :else  (name (key (apply max-key val markmap))))]
    {:cellvalency maxmark :genex genex}))

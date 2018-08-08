(ns bivex.cell
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules]))

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

(defn apply-rule
  "apply the selected rule and update the chromtape"
  [chrom_in]
  (let [prevnuc_idx (chromatin/find-idx-with-head (:chromtape chrom_in))
        nextnuc_idx (chromatin/nucleo-idx-next-head (:chromtape chrom_in))
        prevnuc (chromatin/find-nucleosome-with-head (:chromtape chrom_in))
        rule (rules/select-rule prevnuc (:rules chrom_in))
        prevnuc_new [prevnuc_idx (change-chrom rule prevnuc)] 
        nextnuc_new [nextnuc_idx (move-head (:chromtape chrom_in) nextnuc_idx)]
        new_chromtape (concat
                       (map #(nth (:chromtape chrom_in) %)
                            (get-the-rest-idx prevnuc_idx nextnuc_idx (:chromtape chrom_in)))
                       (vector nextnuc_new)
                       (vector prevnuc_new))
        new_rule (rules/update-rules (:orules chrom_in) nextnuc_new prevnuc_new)]
;    (println nextnuc_new)
    {:k4mono (:k4mono chrom_in)
     :k27mono (:k27mono chrom_in)
     :biv (:biv chrom_in)
     :genex (:genex chrom_in)
     :chromtape (sort new_chromtape)
     :rules new_rule
     :orules (:orules chrom_in)}
    ))

(defn check-valency
  "check valency of the given chromtape"
  [new_chromtape]
  (let [k4total (vec (map #(:k4 (second %)) new_chromtape))
        k27total (vec (map #(:k27 (second %)) new_chromtape))
        biv (map #(bit-and %1 %2) k4total k27total)
        xor (map #(bit-xor %1 %2) k4total k27total)
        k4mono (map #(bit-and %1 %2) k4total xor)
        k27mono (map #(bit-and %1 %2) k27total xor)]
  {:k4mono (vector (apply + k4mono)) :k27mono (vector (apply + k27mono)) :biv (vector (apply + biv))}))

(defn gene-on?
     "calls gene expression outcome"
     [k4mono k27mono biv]
     (cond (> k4mono (+ k27mono 1)          ;(max k27mono biv)
              ) [1]
           :else [0])
     )

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

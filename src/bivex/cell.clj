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
  [no_idx chromtape_idx]
  (remove #(= no_idx %) chromtape_idx))

(defn get-the-rest-idx
  "get the index of the rest of the unchanging nucleosome"
  [prevnuc_idx nextnuc_idx nuchead_idx chromtape]
  (let [sub-chrom (cond (nil? prevnuc_idx) (range (count chromtape))
                        :else (drop-nuc prevnuc_idx (range (count chromtape))))]
    (drop-nuc nuchead_idx (drop-nuc nextnuc_idx sub-chrom))
    ))

(defn move-head
  "assign new head position"
  [chromtape nextnuc_idx]
  (assoc (second (nth chromtape nextnuc_idx)) :head 1))

(defn remove-head
  "remove head for a given nucleosome"
  [nuc]
  (assoc nuc :head 0))

;; issue : when nuchead is 

(defn update_chromtape
  [nuc_c_new nuc_h_new nuc_n_new chromtape]
  (let [rest_n (map #(nth chromtape %)
                    (get-the-rest-idx (first nuc_c_new) (first nuc_n_new) (first nuc_h_new) chromtape))
        temp_chromtape (concat rest_n (vector nuc_c_new) (vector nuc_n_new))]
    (cond (some #(= (first nuc_h_new) (first %)) temp_chromtape) temp_chromtape
          :else (concat temp_chromtape (vector nuc_h_new)))
    ))


(defn apply-rule
  "apply the selected rule and update the chromtape"
  [chrom_in]
  (let [nuc_h (chromatin/find-nucleosome-with-head (:chromtape chrom_in))
        nuc_n_idx (chromatin/nucleo-idx-next-head (:chromtape chrom_in) (first nuc_h))
        rule (rules/select-rule (second nuc_h) (:rules chrom_in))
        nuc_n_new [nuc_n_idx (move-head (:chromtape chrom_in) nuc_n_idx)]
        nuc_c (cond (= (:action rule) "turnover") (chromatin/turnover-match rule (first nuc_h) nuc_n_idx (:chromtape chrom_in))
                    :else nuc_h)
        nuc_c_new (cond (nil? nuc_c) nil :else [(first nuc_c) (change-chrom rule (second nuc_c))])
        nuc_h_new [(first nuc_h) (remove-head (second nuc_h))]
        new_chromtape (update_chromtape nuc_c_new nuc_h_new nuc_n_new (:chromtape chrom_in))
        new_rule (rules/update-rules (:orules chrom_in) nuc_n_new (cond (not= (first nuc_c_new) (first nuc_h_new)) nuc_h_new :else nuc_c_new))
        ]
    (println rule)
  {:k4mono (:k4mono chrom_in)
   :k27mono (:k27mono chrom_in)
   :biv (:biv chrom_in)
   :genex (:genex chrom_in)
   :chromtape (sort new_chromtape)
   :rules (vec new_rule) 
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
     "calls gene expression outcome and multiplies by the transcriptional rate"
     [k4mono k27mono biv trate]
;  (println trate)
  (cond (> k4mono (+ k27mono 1)          ;(max k27mono biv)
           ) [trate]
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

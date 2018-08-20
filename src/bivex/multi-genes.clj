(defn run-multi-genes
  [idx]
  (println idx)
  (spit (clojure.string/join ["output/K27mono-" idx ".txt"])
        (vec (eval/run-bulk (eval/generate_chrom_in @rules/default-rules-file @chromatin/chromatin-file (+ (rand-int 9) 1)) 1000 500 499)))
  )

(pmap run-multi-genes (range 100))

(map #(apply + (read-string (slurp (clojure.string/join ["output/K27mono-" % ".txt"])))) (range 30))

(map #(apply + (read-string (slurp (clojure.string/join ["output/K4mono-" % ".txt"])))) (range 30))

 (jutsu.core/graph! "K27mono distribution" [{:y (map #(apply + (read-string (slurp (clojure.string/join ["output/K27mono-" % ".txt"])))) (range 60)) :type "box"}])

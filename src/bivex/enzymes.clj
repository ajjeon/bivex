(ns bivex.enzymes)


(defn k4methyltransferase
  [nucleo]
  (cond (= (:k4 nucleo) 0) (assoc nucleo :k4 1)
        :else (nucleo))
  )

(defn k27methyltransferase
  [nucleo]
  (cond (= (:k27 nucleo) 0) (assoc nucleo :k27 1)
        :else (nucleo))
  )

(defn k4demethylase
  [nucleo]
  (cond (= (:k4 nucleo) 1) (assoc nucleo :k4 0)
        :else (nucleo))
  )

(defn k4demethylase
  [nucleo]
  (cond (= (:k27 nucleo) 1) (assoc nucleo :k27 0)
        :else ())
  )

(ns bivex.core
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.eval :as eval])
  (:require [jutsu.core :as j])
  (:gen-class))

(defn -main []
  (j/start-jutsu!)
  (println "::: START SINGLE-CELL SIMULATION :::") 
  (eval/run-one chromatin/chromtape 100)
  (println "::: DONE SINGLE-CELL SIMULATION :::")
  (println "::: START BULK SIMULATION :::") 
  (eval/run-bulk chromatin/chromtape 100 1000)
  (println "::: DONE BULK SIMULATION :::") 
)


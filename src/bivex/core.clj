(ns bivex.core
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.eval :as eval])
  (:require [jutsu.core :as j])
  (:require [bivex.rules :as rules])
  (:gen-class))

(def chrom_in {:k4mono [0] :k27mono [0] :biv [0] :genex [0] :chromtape chromatin/chromtape :rules rules/rules}) 

(defn -main []
  (j/start-jutsu!)
  (Thread/sleep 500)
  (println "Opening the browser for graphs...") 
  (println "::: START SINGLE-CELL SIMULATION :::") 
  (eval/run-one chrom_in 100)
  (println "::: DONE SINGLE-CELL SIMULATION :::")
  (println "::: START BULK SIMULATION :::") 
  (eval/run-bulk chrom_in 100 1000)
  (println "::: DONE BULK SIMULATION :::") 
)


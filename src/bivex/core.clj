(ns bivex.core
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.eval :as eval])
  (:require [bivex.cell :as cell])
  (:require [jutsu.core :as j])
  (:require [bivex.rules :as rules])
  (:gen-class))

(def chrom_in
  (let [k4mono (:k4mono (cell/check-valency chromatin/chromtape))
        k27mono (:k27mono (cell/check-valency chromatin/chromtape))
        biv (:biv (cell/check-valency chromatin/chromtape))
        ]
    (print k4mono k27mono biv)
    {:k4mono k4mono
     :k27mono k27mono
     :biv biv
     :genex (cell/gene-on? (first k4mono) (first k27mono) (first biv))
     :chromtape chromatin/chromtape :rules rules/rules}
    )) 

(defn -main [ & args]
  (j/start-jutsu!)
  (Thread/sleep 500)
  (println "Opening the browser for graphs...")
  (println "Reading in rules file")
  (rules/read-in-rules "resources/rules.csv")
  (println "::: START SINGLE-CELL SIMULATION :::") 
  (eval/run-one chrom_in 1000)
  (println "::: DONE SINGLE-CELL SIMULATION :::")
  (println "::: START BULK SIMULATION :::") 
  (eval/run-bulk chrom_in 100 1000)
  (println "::: DONE BULK SIMULATION :::") 
)


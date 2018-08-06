(ns bivex.core
  (:require [clojure.tools.cli :as tools])
  (:require [bivex.eval :as eval])
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.cell :as cell])
  (:require [jutsu.core :as j])
  (:require [bivex.rules :as rules])
  (:require [bivex.files :as files])
  (:gen-class))


(defn generate_chrom_in [rfile cfile]
  (let [chromtape (files/read-in-chromatin cfile)
        k4mono (:k4mono (cell/check-valency chromtape))
        k27mono (:k27mono (cell/check-valency chromtape))
        biv (:biv (cell/check-valency chromtape))
        rules (files/read-in-file rfile)]
;    (print k4mono k27mono biv)
    {:k4mono k4mono
     :k27mono k27mono
     :biv biv
     :genex (cell/gene-on? (first k4mono) (first k27mono) (first biv))
     :chromtape chromtape :rules rules}
    ))

(defn -main
  ""
  [& args]
  (let [[options extra-args banner]
        (tools/cli args
             ["-r" "--rules" "rules"]
             ["-c" "--chromatin" "chromatin"]
             ["-h" "--help" "Show help" :default false :flag true]
             ["-v" "--verbose" "debug mode" :default false :flag true])]

    (reset! rules/default-rules-file (:rules options))
    (reset! chromatin/chromatin-file (:chromatin options))

    (println @rules/default-rules-file)
    (println @chromatin/chromatin-file)
    
    (when (:help options)
      (println banner)
      (System/exit 0))

    (j/start-jutsu!)
    (Thread/sleep 500)
    (println "Opening the browser for graphs...")
    (println "::: START SINGLE-CELL SIMULATION :::") 
    (eval/run-one (generate_chrom_in @rules/default-rules-file @chromatin/chromatin-file) 100)
    (println "::: DONE SINGLE-CELL SIMULATION :::")
    (println "::: START BULK SIMULATION :::") 
    (eval/run-bulk (generate_chrom_in @rules/default-rules-file @chromatin/chromatin-file) 100 100)
    (println "::: DONE BULK SIMULATION :::")
    ))




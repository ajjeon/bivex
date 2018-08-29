(ns bivex.core
  (:require [clojure.tools.cli :as tools])
  (:require [bivex.eval :as eval])
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.cell :as cell])
  (:require [jutsu.core :as j])
  (:require [bivex.rules :as rules])
  (:require [bivex.files :as files])
  (:gen-class))

(def chrom_in (eval/generate_chrom_in @rules/default-rules-file @chromatin/chromatin-file (+ (rand-int 9) 1))) 

(defn -main
  "" ;lein run -b "resources/chromtape.csv" -a "resources/rules.csv" -c 10 -t 1000 -r "resources/new-rules.csv" -n 10
  [& args]
  (let [[options extra-args banner]
        (tools/cli args
             ["-a" "--rules" "rules"]
             ["-b" "--chromatin" "chromatin"]
             ["-c" "--ncells" "ncells"]
             ["-t" "--niters" "niters"]
             ["-r" "--newrules" "newrules"]
             ["-n" "--aftern" "aftern"] ; invoke new rulesets after certain iters
             ["-h" "--help" "Show help" :default false :flag true]
             ["-v" "--verbose" "debug mode" :default false :flag true]
             ["-p" "--plot" "graphs" :default false :flag true])
        beforeiter (read-string (:aftern options)) 
        afteriter (- (read-string (:niters options))
                     (read-string (:aftern options)))]

    (println options)

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

    (when (:plot options)
      (reset! bivex.eval/plot? true))

    (eval/run-one-with-change
     (eval/generate_chrom_in @rules/default-rules-file @chromatin/chromatin-file (+ (rand-int 9) 1))
     beforeiter afteriter)

    (println "::: DONE SINGLE-CELL SIMULATION :::")
    (println "::: START BULK SIMULATION :::")

    (eval/run-bulk
     (read-string (:ncells options))
     beforeiter afteriter)
    
      )
  (println "::: DONE BULK SIMULATION :::")
  (System/exit 0)
    )


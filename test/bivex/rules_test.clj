




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; (defn discourage-biv
;;   [rules givenrules givenm]
;;   (let [changem (cond (= givenm "k4") "k27" :else "k4")
;;         drule (into {} (filter #(and (= (:action %) "methyltransferase")
;;                                       (= (:class %) changem)) givenrules)) 
;;         srule (map #(into {} %) (filter #(or (not= (:action %) "methyltransferase")
;;                                               (not= (:class %) changem)) rules))
;;         new_drule (assoc (assoc drule :affinity 0.5) :abundance 0.5)]
;;     (concat [new_drule] srule)))

;; (defn update-rates
;;   "change feedback ratio"
;;   [drule srule type & {:keys [mfactor]}]
;;   (let [ffactor (cond (nil? mfactor) (cond (= type "methyltransferase") 2 :else 5) :else mfactor) ]
;;     (concat [(assoc (assoc drule
;;                            :affinity (* (:affinity drule) ffactor))
;;                     :abundance (* (:abundance drule) ffactor))] srule)))

;; (defn rule-recruitment
;;   "mimics the recruitment by an existing methyl mark or an empty mark" ;;;;TODO
;;   [givenrules changem type & {:keys [mfactor]}]
;;   (let [drule (into {} (filter #(and (= (:action %) type)
;;                                       (= (:class %) changem)) givenrules)) 
;;         srule (map #(into {} %) (filter #(or (not= (:action %) type)
;;                                              (not= (:class %) changem)) givenrules))]
;; ;    (println type changem)
;; ;    (println givenrules)
;;     (cond (empty? drule) givenrules
;;           :else (update-rates drule srule type :mfactor mfactor))))

;; (defn update-rules-discourage-biv
;;   "based on the new nucleosome, discourage oppositng methyltransferase"
;;   [rules givenrules nextnuc_new]
;;   (let [x (map #(get (second nextnuc_new) %) [:k4 :k27])
;;         xtest (= (apply + x) 1)]
;;     (cond xtest (discourage-biv rules givenrules (name (get-key (+ (.indexOf x 1) 1))))
;;           :else rules)))

;; (defn update-rules-recruitment-sub
;;   "adjust rates accordingly to form loop"
;;   [mark prevnuc_new rules]
;;   (cond (zero? ((keyword mark) (second prevnuc_new)))
;;         (rule-recruitment rules mark "demethylase")
;;         :else (rule-recruitment rules mark "methyltransferase")))

;; (defn update-rules-recruitment
;;   "based on the previous nucleosome, encourage recruitment by existing mark"
;;   [rules prevnuc_new]
;;   (let [updatedrules (update-rules-recruitment-sub "k4" prevnuc_new rules)]
;;     (update-rules-recruitment-sub "k27" prevnuc_new updatedrules)))

;; (defn update-rules-given-locus
;;   "locus-specific recruitment of rules"
;;   [rules nuc]
;;   (cond (= (:locus (second nuc)) 1) (rule-recruitment rules "k4" "methyltransferase" :mfactor 5)
;;         (= (:locus (second nuc)) 2) (rule-recruitment rules "k27" "methyltransferase" :mfactor 5)
;;         :else rules))

;; (defn update-rules
;;   [new_chromtape orules nuc_h]
;;   (let [nextnuc_new (chromatin/find-nucleosome-with-head new_chromtape)
;;         prevnuc_new (nth new_chromtape (first nuc_h))
;;         urules (update-rules-recruitment orules prevnuc_new)
;;         vrules (update-rules-discourage-biv orules urules nextnuc_new)] ;; after every iteration, same default rules get read in
;;     vrules))

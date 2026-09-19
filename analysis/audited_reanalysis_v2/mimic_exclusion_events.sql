SELECT c.stay_id AS id,extract(epoch FROM(c.charttime-i.intime))/3600 AS hr,c.itemid,c.value,c.valuenum
      FROM mimiciv_icu.chartevents c JOIN mimiciv_icu.icustays i USING(stay_id)
      WHERE c.itemid IN(223758,228687,229784,229270,229277,229280)
      AND c.charttime<=i.intime+interval '8 days'
SELECT v.stay_id AS id,
    extract(epoch from (v.starttime-i.intime))/3600 AS start_hr,
    extract(epoch from (v.endtime-i.intime))/3600 AS end_hr,ventilation_status
    FROM mimiciv_derived.ventilation v JOIN mimiciv_icu.icustays i USING(stay_id)
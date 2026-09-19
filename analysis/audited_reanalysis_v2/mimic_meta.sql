SELECT i.stay_id AS id,i.subject_id AS pid,i.hadm_id,
    p.anchor_age+extract(year from a.admittime)-p.anchor_year AS age,p.gender AS sex,
    p.anchor_year_group AS calendar_group,i.first_careunit AS unit,
    extract(epoch from (a.deathtime-i.intime))/3600 AS death_hr,
    extract(epoch from (a.dischtime-i.intime))/3600 AS hosp_end_hr,
    extract(epoch from (i.outtime-i.intime))/3600 AS icu_end_hr,
    a.hospital_expire_flag AS hosp_death
    FROM mimiciv_icu.icustays i JOIN mimiciv_hosp.admissions a USING(hadm_id)
    JOIN mimiciv_hosp.patients p ON p.subject_id=i.subject_id
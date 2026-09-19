SELECT patientunitstayid AS id,uniquepid AS pid,patienthealthsystemstayid AS hadm_id,
    CASE WHEN age ~ '^[0-9]+$' THEN age::numeric WHEN age='> 89' THEN 90 WHEN age='< 18' THEN 17 ELSE NULL END AS age,
    gender AS sex,hospitalid AS site,hospitaldischargeyear AS calendar_group,unittype AS unit,
    hospitaldischargeoffset/60.0 AS hosp_end_hr,unitdischargeoffset/60.0 AS icu_end_hr,
    CASE WHEN hospitaldischargestatus='Expired' THEN 1 WHEN hospitaldischargestatus='Alive' THEN 0 ELSE NULL END AS hosp_death,
    CASE WHEN hospitaldischargestatus='Expired' THEN hospitaldischargeoffset/60.0 ELSE NULL END AS death_hr
    FROM patient
SELECT patientunitstayid AS id,respcarestatusoffset/60.0 AS status_hr,airwaytype,
    ventstartoffset/60.0 AS start_hr,ventendoffset/60.0 AS end_hr FROM respiratorycare
    WHERE airwaytype IS NOT NULL AND airwaytype<>''
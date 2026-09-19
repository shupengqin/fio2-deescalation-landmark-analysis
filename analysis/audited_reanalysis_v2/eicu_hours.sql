WITH fr AS (
      SELECT patientunitstayid AS id,respchartoffset,respchartentryoffset,
      CASE WHEN respchartvalue ~ '^ *[0-9]+([.][0-9]+)? *%? *$'
      THEN replace(trim(respchartvalue),'%','')::numeric ELSE NULL END AS val
      FROM respiratorycharting WHERE respchartvaluelabel IN('FiO2','FIO2 (%)','Set Fraction of Inspired Oxygen (FIO2)')
      AND respchartoffset>=0 AND respchartoffset<11520
    ), f AS (
      SELECT id,floor(respchartoffset/60.0)::int AS hr,
      (array_agg(CASE WHEN val BETWEEN .21 AND 1 THEN val*100 ELSE val END ORDER BY respchartoffset DESC,respchartentryoffset DESC)
      FILTER(WHERE val BETWEEN 21 AND 100 OR val BETWEEN .21 AND 1))[1] AS fio2
      FROM fr GROUP BY id,hr
    ), s AS (
      SELECT patientunitstayid AS id,floor(observationoffset/60.0)::int AS hr,
      (array_agg(sao2 ORDER BY observationoffset DESC) FILTER(WHERE sao2 BETWEEN 1 AND 100))[1] AS sat,
      count(sao2) FILTER(WHERE sao2 BETWEEN 1 AND 100) AS sat_n,
      min(sao2) FILTER(WHERE sao2 BETWEEN 1 AND 100) AS sat_min,
      (array_agg(systemicmean ORDER BY observationoffset DESC) FILTER(WHERE systemicmean BETWEEN 10 AND 250))[1] AS mbp,
      (array_agg(heartrate ORDER BY observationoffset DESC) FILTER(WHERE heartrate BETWEEN 10 AND 300))[1] AS heart_rate,
      (array_agg(respiration ORDER BY observationoffset DESC) FILTER(WHERE respiration BETWEEN 1 AND 100))[1] AS rr_vital
      FROM vitalperiodic WHERE observationoffset>=0 AND observationoffset<11520
      AND patientunitstayid IN(SELECT DISTINCT id FROM fr)
      GROUP BY patientunitstayid,hr
    ) SELECT coalesce(f.id,s.id) AS id,coalesce(f.hr,s.hr) AS hr,fio2,sat,sat_n,sat_min,mbp,heart_rate,rr_vital
    FROM f FULL JOIN s USING(id,hr)
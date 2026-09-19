WITH f AS (
      SELECT v.stay_id AS id,floor(extract(epoch from(v.charttime-i.intime))/3600)::int AS hr,
      (array_agg(fio2 ORDER BY charttime DESC) FILTER(WHERE fio2 BETWEEN 21 AND 100))[1] AS fio2,
      (array_agg(peep ORDER BY charttime DESC) FILTER(WHERE peep BETWEEN 0 AND 40))[1] AS peep,
      (array_agg(respiratory_rate_total ORDER BY charttime DESC) FILTER(WHERE respiratory_rate_total BETWEEN 1 AND 100))[1] AS rr_vent
      FROM mimiciv_derived.ventilator_setting v JOIN mimiciv_icu.icustays i USING(stay_id)
      WHERE v.charttime>=i.intime AND v.charttime<i.outtime AND v.charttime<i.intime+interval '8 days'
      GROUP BY v.stay_id,hr
    ), s AS (
      SELECT v.stay_id AS id,floor(extract(epoch from(v.charttime-i.intime))/3600)::int AS hr,
      (array_agg(spo2 ORDER BY charttime DESC) FILTER(WHERE spo2 BETWEEN 1 AND 100))[1] AS sat,
      count(spo2) FILTER(WHERE spo2 BETWEEN 1 AND 100) AS sat_n,
      min(spo2) FILTER(WHERE spo2 BETWEEN 1 AND 100) AS sat_min,
      (array_agg(coalesce(mbp,mbp_ni) ORDER BY charttime DESC) FILTER(WHERE coalesce(mbp,mbp_ni) BETWEEN 10 AND 250))[1] AS mbp,
      (array_agg(heart_rate ORDER BY charttime DESC) FILTER(WHERE heart_rate BETWEEN 10 AND 300))[1] AS heart_rate,
      (array_agg(resp_rate ORDER BY charttime DESC) FILTER(WHERE resp_rate BETWEEN 1 AND 100))[1] AS rr_vital
      FROM mimiciv_derived.vitalsign v JOIN mimiciv_icu.icustays i USING(stay_id)
      WHERE v.charttime>=i.intime AND v.charttime<i.outtime AND v.charttime<i.intime+interval '8 days'
      GROUP BY v.stay_id,hr
    ) SELECT coalesce(f.id,s.id) AS id,coalesce(f.hr,s.hr) AS hr,fio2,peep,rr_vent,sat,sat_n,sat_min,mbp,heart_rate,rr_vital
    FROM f FULL JOIN s USING(id,hr)
SELECT patientunitstayid AS id,cplitemoffset/60.0 AS hr,cplitemvalue AS value
      FROM public.careplangeneral WHERE cplgroup='Care Limitation' AND cplitemoffset<=8*1440
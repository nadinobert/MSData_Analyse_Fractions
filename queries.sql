-- warum funktioniert die kaskade beim löschen nicht mehr? wenn in analysis oder result gelöscht wird, wird peptide und protein nicht gelöscht

--DELETE FROM analysis
--WHERE id = 107;

-- Get partner proteins RdhA und RdhB
-- set analysis_id of interest
CREATE TEMP VIEW analysis49
AS
   SELECT DISTINCT p.accession, r.sample, p.description
FROM proteins p
    INNER JOIN result r on r.id = p.result_id
    INNER JOIN analysis a on a.id = r.analysis_id
        WHERE analysis_id = 49 AND
              (p.description LIKE '%rdhB%'
                   OR p.description LIKE '%rdhA%');

-- Find RdhB partner for detected RdhA proteins in temp analysis
SELECT sample, p.RdhA_accession, p.RdhB_accession
FROM analysis48 a
INNER JOIN rdhAB_Partner p ON RdhA_accession = a.accession
        WHERE p.RdhB_accession || '&' || a.sample IN (SELECT accession || '&' || sample FROM analysis49)
;

-- Find all distinct RdhA and RdhB accessions in Analysis 48 (a)
SELECT DISTINCT proteins.accession, r.sample, proteins.abundance, proteins.description, proteins.numPeptides
FROM proteins
    INNER JOIN result r on r.id = proteins.result_id
    INNER JOIN analysis a on a.id = r.analysis_id
    WHERE analysis_id = 99 AND
          (proteins.description LIKE '%rdh%')
;

-- count the amount of proteins or peptides that were identified in each sample
SELECT
    r.sample,
    COUNT(DISTINCT proteins.accession) AS protein_count
FROM proteins
    INNER JOIN result r ON r.id = proteins.result_id
    INNER JOIN analysis a ON a.id = r.analysis_id
WHERE a.id = 102
    AND proteins.accession LIKE '%cbdbA%'
GROUP BY r.sample;

-- check if entry contains specific expression
SELECT *
FROM proteins p
WHERE p.accession IN (SELECT RdhB_accession FROM rdhAB_Partner);

-- all proteins from a certain analysis and a certain result wo contaminants
SELECT *
FROM proteins p
    INNER JOIN result r on r.id = p.result_id
    INNER JOIN analysis a on r.analysis_id = a.id
    WHERE
    analysis_id = 109 AND
          p.accession LIKE '%cbdb%' AND
          r.id = 2494
          --r.sample = 'frac11+12' AND
   --(description LIKE '%rdhA%' --AND accession = 'cbdbA0238')
    --OR description LIKE '%rdhB%'
    --OR description LIKE '%OmeA%'
    --OR description LIKE '%OmeB%'
    --OR description LIKE '%hupL%'
    --OR description LIKE '%hupS%'
    --OR description LIKE '%hupX%')
;

-- Top x proteins
-- wo contaminants
-- chose order by t.x (x = desired parameter for ranking)
-- Select only reasonable detected proteins in sample -> top x with highest coverage or abundance, number of unique peptides > x, abundance > x
-- rank function enables to define subgroups (partition) and order subgroups by defined parameter
SELECT * FROM (
    SELECT *, rank() over (
        partition by t.result_id
        order by t.abundance desc
    ) AS rank FROM (
        SELECT
            r.sample,
            proteins.result_id,
            proteins.accession,
            proteins.description,
            proteins.coverage,
            proteins.numPeptides,
            proteins.abundance,
            proteins.MW,
            proteins.numUniquePeptides
        FROM proteins
        INNER JOIN result r on r.id = proteins.result_id
        INNER JOIN analysis a on a.id = r.analysis_id
        WHERE analysis_id = 109 AND
              accession LIKE '%cbdb%' AND
              proteins.abundance <> '' AND
              proteins.numUniquePeptides >= 2 AND
              proteins.abundance > 100000
    ) AS t
) AS u
WHERE u.rank < 11 --AND
   --(description LIKE '%rdhA%' --AND accession = 'cbdbA0238')
    --OR description LIKE '%rdhB%'
    --OR description LIKE '%OmeA%'
    --OR description LIKE '%OmeB%'
    --OR description LIKE '%hupL%'
    --OR description LIKE '%hupS%'
    --OR description LIKE '%hupX%')
;

--Find all (OHR) protein in a specific analysis
--wo contaminants
SELECT *
FROM peptides
    INNER JOIN result r on peptides.result_id = r.id
    INNER JOIN analysis a on a.id = r.analysis_id
    WHERE analysis_id = 95 AND
          peptides.abundance <> ''
      AND peptides.accession LIKE '%cbdbA%'
          --proteins.markedAs = 'False'--AND
             --(description LIKE '%rdhA%'
        --OR description LIKE '%rdhB%'
        --OR description LIKE '%OmeA%'
        --OR description LIKE '%OmeB%'
        --OR description LIKE '%hupL%'
        --OR description LIKE '%hupS%'
        --OR description LIKE '%hupX%')
ORDER BY peptides.result_id, peptides.abundance DESC
;
--check a set of results for contaminants -> "marked as" '' (nothing) or markedAs True
SELECT *
FROM proteins
    INNER JOIN result r on proteins.result_id = r.id
    INNER JOIN analysis a on a.id = r.analysis_id
    WHERE analysis_id = 85 AND proteins.markedAs = 'False'
ORDER BY proteins.result_id, proteins.abundance DESC
;

-- Suche cbdbB003 peptide in bestimmten ansätzen
SELECT p.confidence, p.sequence, p.modifications, p.numPSMs, p.accession, p.abundance, p.xCorr, r.sample, a.date
FROM peptides p
INNER JOIN result r on p.result_id = r.id
INNER JOIN analysis a on a.id = r.analysis_id
where p.accession = 'cbdbB0003'
    AND (result_id = 77
    OR result_id = 78
    OR result_id = 81
    OR result_id = 82)
;

-- Find approaches where aminoeptidase cbdbA1001, cbdbA1390, cbdbA1394, cbdbA0861 was detected
-- in approaches where membrane fraction was seperated (>49)
SELECT p.confidence, p.sequence, p.accession, p.abundance, r.sample, a.date, r.analysis_id
FROM peptides p
INNER JOIN result r on p.result_id = r.id
INNER JOIN analysis a on a.id = r.analysis_id
where (p.accession = 'cbdbA1001'
      OR p.accession = 'cbdbA1390'
      OR p.accession = 'cbdbA1394'
      OR p.accession = 'cbdbA0861')
    AND r.analysis_id >= 49
    AND p.abundance >= 1000000
;

-- Find samples in one analysis where Keratin was detected
SELECT DISTINCT *
FROM proteins p
INNER JOIN result r on p.result_id = r.id
INNER JOIN analysis a on a.id = r.analysis_id
where (p.description LIKE '%Keratin%'
    OR p.description LIKE '%KERATIN%')
    AND r.analysis_id = 91
;

-- overview of proteins and related peptides in analysis (extendable to ohr-proteins only)
SELECT DISTINCT q.accession, r.sample, q.description, q.abundance, p.sequence, p.abundance, r.RdhB_accession, q.coverage
FROM proteins q
INNER JOIN result r on r.id = q.result_id
INNER JOIN peptides p on q.result_id = p.result_id and q.accession = p.accession
LEFT JOIN rdhAB_Partner r on p.accession = RdhA_accession
INNER JOIN analysis a on r.analysis_id = a.id
WHERE
   a.id = 83 AND
   (description LIKE '%rdhA%'
    OR description LIKE '%rdhB%'
    OR description LIKE '%OmeA%'
    OR description LIKE '%OmeB%'
    OR description LIKE '%hupL%'
    OR description LIKE '%hupS%'
    OR description LIKE '%hupX%')
;

-- Number of unique detected peptides/ proteins per sample respectely (extendable to ohr-proteins only)
select a.date, a.number, r.sample, COUNT(DISTINCT q.accession) AS '#proteine', COUNT(DISTINCT p.sequence) AS '#peptide'
from proteins q
inner join result r on r.id = q.result_id
inner join peptides p on q.result_id = p.result_id and q.accession = p.accession
inner join analysis a on r.analysis_id = a.id

WHERE
a.id = 68
--   AND
--   (description LIKE '%rdhA%'
--   OR description LIKE '%rdhB%'
--   OR description LIKE '%OmeA%'
--   OR description LIKE '%OmeB%'
--   OR description LIKE '%hupL%'
--   OR description LIKE '%hupS%'
--   OR description LIKE '%hupX%')
GROUP BY r.sample
;

-- Create new table for RdhA-RdhB-Partners
--CREATE TABLE rdhABC_Partner
--(
--    RdhA_accession TEXT PRIMARY KEY,
--    RdhB_accession TEXT
--)
--;

-- Insert Partner proteins (RdhA, RdhB) partner looked up at KEGG genome CBDB1
--INSERT INTO rdhAB_Partner (RdhA_accession, RdhB_accession)
--VALUES
       --('cbdbA0080', 'cbdbB0003'),
       --('cbdbA0084', 'cbdbA0085'),
       --('cbdbA0088', 'cbdbB0004'),
       --('cbdbA0096', 'cbdbB0005'),
       --('cbdbA0187', 'cbdbA0188'),
       --('cbdbA0238', 'cbdbA0239'),
       --('cbdbA0243', 'cbdbC0001'),
       --('cbdbA1092', 'cbdbA1094'),
       --('cbdbA1453', 'cbdbA1452'),
       --('cbdbA1455', 'cbdbA1454'),
       --('cbdbA1491', 'cbdbA1490'),
       --('cbdbA1495', 'cbdbB0033'),
       --('cbdbA1503', 'cbdbA1502'),
       --('cbdbA1508', 'cbdbA1507'),
       --('cbdbA1535', 'cbdbA1536'),
       --('cbdbA1539', 'cbdbB0040'),
       --('cbdbA1542', 'cbdbA1541'),
       --('cbdbA1546', 'cbdbA1545'),
       --('cbdbA1550', 'cbdbA1549'),
       --('cbdbA1560', 'cbdbA1559'),
       --('cbdbA1563', 'cbdbA1562'),
       --('cbdbA1570', 'cbdbA1569'),
       --('cbdbA1575', 'cbdbA1573'),
       --('cbdbA1578', 'cbdbA1577'),
       --('cbdbA1582', 'cbdbA1581'),
       --('cbdbA1588', 'cbdbA1587'),
       --('cbdbA1595', 'cbdbA1594'),
       --('cbdbA1598', 'cbdbA1597'),
       --('cbdbA1618', 'cbdbA1617'),
       --('cbdbA1624', 'cbdbA1623'),
       --('cbdbA1627', 'cbdbA1626'),
       --('cbdbA1638', 'cbdbA1637')
       --;
--SELECT * FROM transmembrane_areas WHERE accession = 'cbdbA0131';

--UPDATE rdhAB_Partner SET RdhB_accession = 'cbdbA1502'
--WHERE RdhA_accession = 'cbdbA1503';

-- show only proteins with CAI-accesion no
SELECT DISTINCT accession
FROM proteins
WHERE accession LIKE '%CAI%'
;

-- update proteins with CAI-accession no. to usual cbdbA... - accessions
-- accession is foreign key
UPDATE proteins SET accession = 'cbdbA0130'
WHERE  accession = 'CAI82388.1'
;

-- select all detected Proteins containing TM areas in specific approach (set analysisID's)!
-- abundance >10E5
-- file name: TMProtein
SELECT DISTINCT p.accession
FROM proteins p
INNER JOIN result r
    ON p.result_id = r.id
INNER JOIN analysis a
    on a.id = r.analysis_id
INNER JOIN transmembrane_areas ta
    ON p.accession = ta.accession
WHERE a.id = 57
  --AND abundance > 100000 AND p.abundance != ''
;

-- select distinct proteins in certain analysis
SELECT DISTINCT p.accession
FROM proteins p
INNER JOIN result r on r.id = p.result_id
INNER JOIN analysis a on a.id = r.analysis_id
WHERE a.id = 94
;

-- selects all detected proteins without TM area
-- --with abundance > 10E6 (set analysis ID))
SELECT DISTINCT p.accession
FROM proteins p
INNER JOIN result r
    ON p.result_id = r.id
INNER JOIN analysis a
    ON a.id = r.analysis_id
LEFT JOIN transmembrane_areas ta ON p.accession = ta.accession
WHERE a.id = 57 AND ta.accession IS NULL
  --AND abundance > 1000000 AND p.abundance != ''
;

-- selects all detected proteins with TM area
-- no contaminants
-- (with abundance > 10E6)
-- set analysis ID
SELECT DISTINCT p.accession
FROM proteins p
INNER JOIN result r
    ON p.result_id = r.id
INNER JOIN analysis a
    ON a.id = r.analysis_id
LEFT JOIN transmembrane_areas ta ON p.accession = ta.accession
WHERE a.id = 84
  AND ta.accession IS NOT NULL
  AND p.accession LIKE '%cbdbA%'
  --AND abundance > 1000000 AND p.abundance != ''
;

-- selects all detected proteins with TM area
-- under the top 10 coverage/abundance/unique peptides ect per analysis
-- rank function enables to define subgroups (partition) and order subgroups by defined parameter
-- order by abundance, coverage or amount of peptides
-- no contaminants
-- $
SELECT DISTINCT * FROM (
    SELECT *, rank() over (
        partition by t.result_id
        order by t.abundance desc
    ) AS rank FROM (
        SELECT
            r.sample,
            proteins.result_id,
            proteins.accession,
            proteins.description,
            proteins.coverage,
            proteins.numPeptides,
            proteins.abundance,
            proteins.MW
        FROM proteins
        INNER JOIN result r on r.id = proteins.result_id
        INNER JOIN analysis a on a.id = r.analysis_id
        LEFT JOIN transmembrane_areas ta ON proteins.accession = ta.accession
WHERE a.id = 92
    AND ta.accession IS NOT NULL
    --AND proteins.accession LIKE '%cbdbA%'
    AND proteins.abundance > 100000 AND proteins.numPeptides > 2
    ) AS t
) AS u
WHERE u.rank < 11
;


-- select peptides with dimethylatet/ acetylated amino acids
SELECT *
FROM peptides t
WHERE t.result_id = 1546
  AND t.modifications LIKE '%Acetyl%'
  --AND t.confidence = 'High'--AND
  --(description LIKE '%rdhA%'
  --OR description LIKE '%rdhB%'
  --OR description LIKE '%OmeA%'
  --OR description LIKE '%OmeB%'
  --OR description LIKE '%hupL%'
  --OR description LIKE '%hupS%'
  --OR description LIKE '%hupX%')

;

-- query.sql
WITH RECURSIVE split(result_id, masterModification, str) AS (
    SELECT result_id, '', wurm||';' FROM (

       SELECT x.result_id, x.masterModifications, masterModifications as wurm FROM
(SELECT t.result_id, t.masterModifications, instr(t.masterModifications, '[') AS a, instr(t.masterModifications, ']') AS b
FROM peptides t
WHERE t.result_id = 1447
  AND t.modifications LIKE '%Acetyl%'
  AND t.confidence = 'High'--AND
  --(description LIKE '%rdhA%'
  --OR description LIKE '%rdhB%'
  --OR description LIKE '%OmeA%'
  --OR description LIKE '%OmeB%'
  --OR description LIKE '%hupL%'
  --OR description LIKE '%hupS%'
  --OR description LIKE '%hupX%')
    ) x
                                                        )
    UNION ALL SELECT
    result_id,
    substr(str, 0, instr(str, ';')),
    substr(str, instr(str, ';')+1)
    FROM split WHERE str!=''
)
SELECT result_id, masterModification
FROM split WHERE result_id = 1446;


-- shows all Proteins and peptides when TM areas were detected in specific analysis ID's (set analysis ID's)
-- file name: TMAreas
select * from
transmembrane_areas t
inner join
(select x.sequence, x.result_id, x.groupId, x.accession, substr(x.pos, x.a+1, x.b-x.a-1) AS start, substr(x.pos, x.b+1, x.c-x.b-1) AS end from
(select p.sequence, p.result_id, p.groupId, p.accession, p.position AS pos, instr(p.position, '[') AS a, instr(p.position, '-') AS b,
        instr(p.position, ']') AS c FROM peptides p) x) y
ON t.accession = y.accession

inner join
    result r on r.id = y.result_id
inner join
    analysis a2 on r.analysis_id = a2.id
WHERE ((y.start >= t.start AND y.start <= t.end) OR (y.end >= t.start AND y.end <= t.end))
  AND
      (a2.id = 104)
;


-- extend table with detected TM areas
-- shows all proteins and peptides with detected TM areas in specific analysis ID's (set analysis ID's)
-- + shows if the detected TM area was targeted by modification
-- file to save for marie name: TMAreas2
SELECT DISTINCT t.accession, t.start, t.end,
                y.result_id, y.sequence, y.start, y.end,
                p2.result_id, p2.modifications, p2.position,
                r.id, r.sample
FROM
transmembrane_areas t
INNER JOIN
(SELECT x.sequence, x.result_id, x.groupId, x.accession, substr(x.pos, x.a+1, x.b-x.a-1) AS start, substr(x.pos, x.b+1, x.c-x.b-1) AS end FROM
(SELECT p.sequence, p.result_id, p.groupId, p.accession, p.position AS pos, instr(p.position, '[') AS a, instr(p.position, '-') AS b,
        instr(p.position, ']') AS c FROM peptides p) x) y
ON t.accession = y.accession

INNER JOIN
    result r ON r.id = y.result_id
INNER JOIN
    peptides p2 ON p2.result_id = y.result_id
    AND
           p2.sequence = y.sequence
    AND
           p2.accession = t.accession
INNER JOIN
    analysis a2 ON r.analysis_id = a2.id
WHERE ((y.start >= t.start AND y.start <= t.end)
           OR
       (y.end >= t.start AND y.end <= t.end))
  AND
      (a2.id = 34 OR a2.id = 35 OR a2.id = 36)
  AND
      p2.modifications <> ''
;

-- check all ever detected proteins for proteases
SELECT DISTINCT p.description, r.id, r.sample, a.date, a.comment
FROM proteins p
INNER JOIN result r on p.result_id = r.id
INNER JOIN analysis a on a.id = r.analysis_id
WHERE p.description LIKE '%protease%'
;

-- shows if peptides from a certain protein were detected over all trials
-- file name: TMAreas
select * from
peptides p
WHERE p.accession = 'cbdbA1690'
;

-- try out chat gpt code for Dissertation plotting:
SELECT * FROM (
    SELECT *, RANK() OVER (
        PARTITION BY t.result_id
        ORDER BY t.abundance DESC
    ) AS rank FROM (
        SELECT
            r.sample,
            proteins.result_id,
            proteins.accession,
            proteins.description,
            proteins.coverage,
            proteins.numPeptides,
            proteins.abundance,
            proteins.MW,
            proteins.numUniquePeptides
        FROM proteins
        INNER JOIN result r ON r.id = proteins.result_id
        INNER JOIN analysis a ON a.id = r.analysis_id
        WHERE r.analysis_id = 92
            AND proteins.accession LIKE 'cbdbA%'
            AND proteins.abundance <> ''
            AND proteins.numUniquePeptides >= 2
            AND proteins.abundance > 100000
    ) AS t
) AS u
WHERE u.rank < 50
    AND (
        u.description LIKE '%rdhA%' -- AND u.accession = 'cbdbA0238'
        OR u.description LIKE '%rdhB%'
        OR u.description LIKE '%OmeA%'
        OR u.description LIKE '%OmeB%'
        OR u.description LIKE '%hupL%'
        OR u.description LIKE '%hupS%'
        OR u.description LIKE '%hupX%'
    )
UNION ALL

SELECT
    r.sample,
    proteins.result_id,
    proteins.accession,
    proteins.description,
    proteins.coverage,
    proteins.numPeptides,
    proteins.abundance,
    proteins.MW,
    proteins.numUniquePeptides,
    NULL as rank
FROM proteins
INNER JOIN result r ON r.id = proteins.result_id
INNER JOIN analysis a ON a.id = r.analysis_id
WHERE r.analysis_id = 92
    AND proteins.accession LIKE 'cbdbA%'
    AND proteins.abundance <> ''
    AND (proteins.description LIKE '%rdhB%'
        OR proteins.description LIKE '%omeB%')
;
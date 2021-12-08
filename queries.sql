
SELECT *
FROM peptides t
    INNER JOIN result r on t.result_id = r.id
        WHERE t.accession = "cbdbA0193"
;

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
    WHERE analysis_id = 48 AND
          (proteins.description LIKE '%rdhA%')
;

-- check if entry contains specific expression
SELECT *
FROM proteins p
WHERE p.accession IN (SELECT RdhB_accession FROM rdhAB_Partner);

-- analyse SDS bands Top 5
-- Select only reasonable detected proteins in gel bands -> top 5 with highest coverage, top 5 with highest abundance (<6)
-- rank function enables to define subgroups (partition) and order subgroups by defined parameter
SELECT * FROM (
    SELECT *, rank() over (
        partition by t.result_id
        order by t.coverage desc
    ) AS rank FROM (
        SELECT
            proteins.result_id,
            proteins.accession,
            proteins.description,
            proteins.coverage,
            proteins.numPeptides,
            proteins.abundance,
            proteins.MW,
            r.sample
        FROM proteins
        INNER JOIN result r on r.id = proteins.result_id
        INNER JOIN analysis a on a.id = r.analysis_id
        WHERE analysis_id = 49 AND proteins.abundance <> ''
    ) AS t
) AS u
WHERE u.rank < 6 AND
   (description LIKE '%rdhA%'
    OR description LIKE '%rdhB%'
    OR description LIKE '%OmeA%'
    OR description LIKE '%OmeB%'
    OR description LIKE '%hupL%'
    OR description LIKE '%hupS%'
    OR description LIKE '%hupX%')
;

SELECT *
FROM proteins
    INNER JOIN result r on proteins.result_id = r.id
    INNER JOIN analysis a on a.id = r.analysis_id
    WHERE analysis_id = 49 AND proteins.abundance <> '' AND r.sample = 18.6
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

-- overview of proteins and related peptides in database (extendable to ohr-proteins only)
SELECT DISTINCT q.accession, r.sample, q.description, q.abundance, p.sequence, p.abundance, r.RdhB_accession, q.coverage
FROM proteins q
INNER JOIN result r on r.id = q.result_id
INNER JOIN peptides p on q.result_id = p.result_id and q.accession = p.accession
LEFT JOIN rdhAB_Partner r on p.accession = RdhA_accession
INNER JOIN analysis a on r.analysis_id = a.id
WHERE
   a.id = 49 AND
   (description LIKE '%rdhA%'
    OR description LIKE '%rdhB%'
    OR description LIKE '%OmeA%'
    OR description LIKE '%OmeB%'
    OR description LIKE '%hupL%'
    OR description LIKE '%hupS%'
    OR description LIKE '%hupX%')
;

-- Amount of unique detected peptides/ proteins per sample (extendable to ohr-proteins only)
select a.date, a.number, r.sample, COUNT(DISTINCT q.accession) AS '#proteine', COUNT(DISTINCT p.sequence) AS '#peptide'
from proteins q
inner join result r on r.id = q.result_id
inner join peptides p on q.result_id = p.result_id and q.accession = p.accession
inner join analysis a on r.analysis_id = a.id

where
   a.date = '2021-11-24' AND
   (description LIKE '%rdhA%'
   OR description LIKE '%rdhB%'
   OR description LIKE '%OmeA%'
   OR description LIKE '%OmeB%'
   OR description LIKE '%hupL%'
   OR description LIKE '%hupS%'
   OR description LIKE '%hupX%')
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

--DELETE FROM analysis
--WHERE id = 46;

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
-- file name: TMProtein
SELECT DISTINCT p.accession, r.sample, a.id, p.abundance, p.confidenceSample
FROM proteins p
INNER JOIN result r
    ON p.result_id = r.id
INNER JOIN analysis a
    on a.id = r.analysis_id
INNER JOIN transmembrane_areas ta
    ON p.accession = ta.accession
WHERE ((a.id = 34 OR a.id = 35 OR a.id = 36))
;

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
      (a2.id = 34 OR a2.id = 35 OR a2.id = 36)
;

-- extend table with detected TM areas
-- shows all proteins and peptides with detected TM areas in specific analysis ID's (set analysis ID's)
-- + shows if the detected TM area was targeted by modification
-- file name: TMAreas2
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

--SELECT * FROM proteins p
--           WHERE p.modifications <> ''
--;
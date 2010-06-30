/* SQL script that modifies the CHADO schema to include tsvector columns for 
   full-text searching by the Bio::DB::Das::Chado adaptor (also requires the 
   modified Chado.pm and Segment.pm modules to work).

   o The script first drops, then creates the searchable tsvector columns 
   representing each target column from the search.  We also populate the
   column from existing data.
   o A trigger function is then added to each table so that any INSERT or 
   UPDATE operations on the target columns repopulate the searchable 
   column.
   o Finally for each table, a new index is created on the searchable column. 
   o Then we recreate the all_feature_names VIEW with a searchable column

   For the trigger functions, we're assuming that we don't modify any existing 
   scripts or modules (other than the Chado adaptor search functions) to
   accommodate our schema change.  This is perhaps slower than modifying all 
   those scripts, but for now it's likely to be more robust, and quicker to
   implement.

   Leighton Pritchard 2010
*/


/* 1) Drop and create searchable target columns */
-- feature table
ALTER TABLE feature DROP COLUMN searchable_name;
ALTER TABLE feature ADD COLUMN searchable_name tsvector;
UPDATE feature SET searchable_name = to_tsvector('pg_catalog.english', COALESCE((name,'') || ' ' || (uniquename,'')));
-- synonym table
ALTER TABLE synonym DROP COLUMN searchable_synonym_sgml;
ALTER TABLE synonym ADD COLUMN searchable_synonym_sgml tsvector;
UPDATE synonym SET searchable_synonym_sgml = to_tsvector('pg_catalog.english', synonym_sgml);
-- dbxref table
ALTER TABLE dbxref DROP COLUMN searchable_accession;
ALTER TABLE dbxref ADD COLUMN searchable_accession tsvector;
UPDATE dbxref SET searchable_accession = to_tsvector('pg_catalog.english', accession);

/* 2) Add trigger function to each table to populate the
      searchable column when a data-modifying operation occurs 
      on the target field 

     This is made much easier by the existence of the tsvector_update_trigger() 
     procedure
*/
-- feature table
DROP TRIGGER IF EXISTS feature_searchable_iu ON feature;
CREATE TRIGGER feature_searchable_iu
  BEFORE INSERT OR UPDATE ON feature
    FOR EACH ROW 
      EXECUTE PROCEDURE 
        tsvector_update_trigger(searchable_name, 'pg_catalog.english', name, uniquename);

-- synonym table
DROP TRIGGER IF EXISTS synonym_searchable_iu ON synonym;
CREATE TRIGGER synonym_searchable_iu
  BEFORE INSERT OR UPDATE ON synonym
    FOR EACH ROW 
      EXECUTE PROCEDURE 
        tsvector_update_trigger(searchable_synonym_sgml, 'pg_catalog.english', synonym_sgml);

-- synonym table
DROP TRIGGER IF EXISTS dbxref_searchable_iu ON dbxref;
CREATE TRIGGER dbxref_searchable_iu
  BEFORE INSERT OR UPDATE ON dbxref
    FOR EACH ROW 
      EXECUTE PROCEDURE 
        tsvector_update_trigger(searchable_accession, 'pg_catalog.english', accession);

/* 3) We need to recreate the all_feature_names VIEW with a searchable 
      column for the name data.  The responses for gmod_materialized_view_tool
      need to be modified.  Presumably, it only makes sense to do this
      with a materialzed view, as without it the performance will no doubt
      be worse.

-- To materialize this view, run gmod_materialized_view_tool.pl -c and
-- answer the questions with these responses:
--
--   all_feature_names
--
--   public.all_feature_names
--
--   y   (yes, replace the existing view)
--
--   (some update frequency, I chose daily)
--
--   feature_id integer,name varchar(255),organism_id integer,searchable_name tsvector
--
--   (the select part of the view below, all on one line)
--
--   feature_id,name
--
--   create index all_feature_names_lower_name on all_feature_names (lower(name))
--
--   y 
*/

CREATE OR REPLACE VIEW all_feature_names (
  feature_id,
  name,
  organism_id,
  searchable_name
) AS
SELECT feature_id, CAST(substring(uniquename FROM 0 FOR 255) AS varchar(255)) AS name, 
       organism_id, to_tsvector('english', CAST(substring(uniquename FROM 0 FOR 255) AS varchar(255))) 
       FROM feature
  UNION
SELECT feature_id, name, organism_id, to_tsvector('english', name) 
       FROM feature 
         WHERE name IS NOT NULL
  UNION
SELECT fs.feature_id, s.name, f.organism_id, to_tsvector('english', s.name) 
       FROM feature_synonym fs, synonym s, feature f
         WHERE fs.synonym_id = s.synonym_id AND fs.feature_id = f.feature_id
  UNION
SELECT fp.feature_id, CAST(substring(fp.value FROM 0 FOR 255) AS varchar(255)) AS name,
       f.organism_id, to_tsvector('english',CAST(substring(fp.value FROM 0 FOR 255) AS varchar(255))) 
       FROM featureprop fp, feature f
         WHERE f.feature_id = fp.feature_id
  UNION
SELECT fd.feature_id, d.accession, f.organism_id,to_tsvector('english',d.accession) 
       FROM feature_dbxref fd, dbxref d,feature f
         WHERE fd.dbxref_id = d.dbxref_id AND fd.feature_id = f.feature_id;


/* 4) Now we index the searchable columns for each modified table.
      The DROP INDEX is a precaution, it will have been dropped along
      with the column it was indexing, in step 1.  That's why we don't 
      check with IF EXISTS, first.
*/
-- feature table
-- DROP INDEX searchable_feature_name_idx;
CREATE INDEX searchable_feature_name_idx ON feature USING gin(searchable_name);
-- synonym table
-- DROP INDEX searchable_synonym_sgml_idx;
CREATE INDEX searchable_synonym_sgml_idx ON synonym USING gin(searchable_synonym_sgml);
-- dbxref table
-- DROP INDEX searchable_dbxref_accession_idx;
CREATE INDEX searchable_dbxref_accession_idx ON dbxref USING gin(searchable_accession);

--all_feature_names materialzed view
CREATE INDEX searchable_all_feature_names_idx ON all_feature_names USING gin(searchable_name);


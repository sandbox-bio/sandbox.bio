#!/bin/bash

pg_dump --clean --if-exists --quote-all-identifiers -h db.bqjvxpdzkembvixymfae.supabase.co -U postgres > dev.sql
pg_dump --clean --if-exists --quote-all-identifiers -h db.jpdymnmaakzeyqyfomcs.supabase.co -U postgres > stg.sql

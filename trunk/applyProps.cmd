REM copy this script in the root directory of the project
REM and run it: ".\applyProps.cmd"

for /r . %%X in (*.py) do (
svn propset svn:eol-style native "%%X"
svn propset svn:keywords "LastChangedBy LastChangedDate Id LastChangedRevision URL" "%%X"
svn propset svn:mime-type text/x-python "%%X"
svn propset svn:executable ON "%%X"
)

for /r . %%X in (*.csv *.rst *.cfg AUTHORS CHANGES LICENSE TODO) do (
svn propset svn:eol-style native "%%X"
svn propset svn:keywords "LastChangedBy LastChangedDate Id LastChangedRevision URL" "%%X"
svn propset svn:mime-type text/plain "%%X"
)



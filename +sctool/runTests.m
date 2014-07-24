function result = runTests

import matlab.unittest.TestSuite;

dirname = which('lapsolve');
[pathstr,fname,fext] = fileparts(dirname);
dirname = fullfile(pathstr,'tests');
suiteClass = TestSuite.fromFolder(dirname);
result = run(suiteClass);

end
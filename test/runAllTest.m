baseFolder = 'Z:/MatlabProjs/LP-v4';
testBaseFolder = strcat(baseFolder, '/test');
suiteFolder = testsuite(...
    {strcat(testBaseFolder, '/lp4'),...
    strcat(testBaseFolder, '/util')});
result = run(suiteFolder);

classdef testStrip < matlab.unittest.TestCase
    
    properties
        map
    end
    
    methods (TestClassSetup)
        function createMap(testCase)
            opt = sctool.scmapopt('trace',0,'tol',1e-12);
            p = polygon([4 2i -2+4i -3 -3-1i 2-2i]);
            testCase.map = stripmap(p,[1 4],opt);
        end
    end
    
    methods (Test)
        
        function testForwardMap(testCase)
            result = testCase.map([1+1i Inf -2+0.5i 0].');
            expected = [
                -3.000000000000000 - 0.312856946533931i
                -3.000000000000000 + 0.000000000000000i
                3.534971699300866 - 0.074755569347609i
                0.000000000000000 + 2.000000000000000i
                ];
            testCase.verifyEqual(result,expected,'abstol',1e-10);
        end
        
        function testInverseMap(testCase)
            val = testCase.map([4 -2+0.5i 1+0.5i]);
            result = evalinv( testCase.map, val );
            expected = [4 -2+0.5i 1+0.5i];
            testCase.verifyEqual(result,expected,'abstol',1e-10);
        end
        
        function testPlot(testCase)
            fig = figure;
            plot(testCase.map,4,3)
            close(fig)
        end
        
    end
    
    
    
end
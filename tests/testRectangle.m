classdef testRectangle < matlab.unittest.TestCase
    
    properties
        map
    end
    
    methods (TestClassSetup)
        function createMap(testCase)
            opt = sctool.scmapopt('trace',0,'tol',1e-12);
            p = polygon([4 2i -2+4i -3 -3-1i 2-2i]);
            testCase.map = rectmap(p,1:4,opt);
        end
    end
    
    methods (Test)
        
        function testForwardMap(testCase)
            result = testCase.map([1.5 1.4+3i -0.6+1i 1].');
            expected = [
                3.641550444027862 - 0.358449555972138i
                -0.005336970451055 + 1.988135598448822i
                -1.643459212104280 + 0.428597577267735i
                1.646072527976422 - 1.929214505595285i
                ];
            testCase.verifyEqual(result,expected,'abstol',1e-7);
        end
        
        function testInverseMap(testCase)
            val = testCase.map([1.5 1.4+3i -0.6+1i 1]);
            result = evalinv( testCase.map, val );
            expected = [1.5 1.4+3i -0.6+1i 1];
            testCase.verifyEqual(result,expected,'abstol',1e-7);
        end
                
        function testPlot(testCase)
            fig = figure;
            plot(testCase.map,4,3)
            close(fig)
        end
        
    end
    
    
    
end
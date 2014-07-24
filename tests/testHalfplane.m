classdef testHalfplane < matlab.unittest.TestCase
    
    properties
        map
    end
   
    methods (TestClassSetup)
        function createMap(testCase)
            opt = sctool.scmapopt('trace',0,'tol',1e-12);
            p = polygon([4 2i -2+4i -3 -3-1i 2-2i]);
            testCase.map = hplmap(p,opt);
        end
    end
    
    methods (Test)
        
        function testForwardMap(testCase)
            result = testCase.map([-1+0.01i 2i 4+0.5i Inf]);
            expected = [
                3.718839996085665 - 0.046084791699413i,...
                1.734612962216089 - 0.777136490010106i,...
                1.179285609480821 - 1.753573737693204i,...
                2.000000000000000 - 2.000000000000000i,...
                ];
            testCase.verifyEqual(result,expected,'abstol',1e-10);
        end
        
        function testInverseMap(testCase)
            val = testCase.map([-1 -1+0.01i 2i 4+0.5i]);
            result = evalinv( testCase.map, val );
            expected = [ -1 -1+0.01i 2i 4+0.5i ];
            testCase.verifyEqual(result,expected,'abstol',1e-10);
        end
        
        function testPlot(testCase)
            fig = figure;
            plot(testCase.map,2,3)
            close(fig)
        end
        
    end
    
    
    
end
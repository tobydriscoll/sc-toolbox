classdef testDisk < matlab.unittest.TestCase
    
    properties
        map
    end
    
    
    %     methods(TestMethodSetup)
    %             function createFigure(testCase)
    %             end
    %     end
    %
    %     methods(TestMethodTeardown)
    %     end
    
    methods (TestClassSetup)
        function createMap(testCase)
            opt = sctool.scmapopt('trace',0,'tol',1e-12);
            p = polygon([4 2i -2+4i -3 -3-1i 2-2i]);
            testCase.map = diskmap(p,opt);
        end
    end
    
    methods (Test)
        
        function testForwardMap(testCase)
            result = testCase.map([0.5+0.5i -0.9 -0.8+0.3i (1+1i)/sqrt(2)]);
            expected = [
                -2.301479291389453 + 0.891455618349974i,...
                -2.959017053517382 - 0.004724964608807i,...
                -2.920229222824237 + 0.110570172682907i,...
                -2.699042997340806 + 1.203828010636846i,...
                ];
            testCase.verifyEqual(result,expected,'abstol',1e-10);
        end
        
        function testInverseMap(testCase)
            val = testCase.map([0.5+0.5i -0.9 -0.8+0.3i (1+1i)/sqrt(2)]);
            result = evalinv( testCase.map, val );
            expected = [0.5+0.5i -0.9 -0.8+0.3i (1+1i)/sqrt(2)];
            testCase.verifyEqual(result,expected,'abstol',1e-10);
        end
        
        function testCenter(testCase)
            f = center(testCase.map,0);
            result = f(0);
            testCase.verifyEqual(result,0,'abstol',1e-11);
        end
        
        function testPlot(testCase)
            fig = figure;
            plot(testCase.map,4,3)
            close(fig)
        end
        
    end
    
    
    
end
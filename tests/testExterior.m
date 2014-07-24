classdef testExterior < matlab.unittest.TestCase
    
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
            testCase.map = extermap(p,opt);
        end
    end
    
    methods (Test)
        
        function testForwardMap(testCase)
            result = testCase.map([0.5+0.5i -0.9 -0.8+0.3i (1+1i)/sqrt(2)]).';
            expected = [
                3.383320944105799 - 2.338017988574543i
                -3.257095120413423 + 0.536117032603197i
                -3.428512520416633 - 0.641446812228358i
                2.571844815094428 - 1.428155184905573i
                ];
            testCase.verifyEqual(result,expected,'abstol',1e-10);
        end
        
        function testInverseMap(testCase)
            val = testCase.map([0.5+0.5i -0.9 -0.8+0.3i (1+1i)/sqrt(2)]);
            result = evalinv( testCase.map, val );
            expected = [0.5+0.5i -0.9 -0.8+0.3i (1+1i)/sqrt(2)];
            testCase.verifyEqual(result,expected,'abstol',1e-10);
        end
        
        function testPlot(testCase)
            fig = figure;
            plot(testCase.map,4,3)
            close(fig)
        end
        
    end
    
    
    
end
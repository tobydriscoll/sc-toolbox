classdef testAnnulus < matlab.unittest.TestCase
    
    properties
        map
    end
    
    methods (TestClassSetup)
        function createMap(testCase)
            p0 = polygon([-2.0335 + 3.0056i  -3.3743 - 0.4358i  -2.3911 - 2.5363i   0.3799 - 3.2737i   2.6816 - 1.3073i   3.2402 + 2.6704i   1.1844 + 3.3855i]);
            p1 = polygon([ -0.5363 + 0.8603i  -0.6034 + 0.2570i   0.2682 - 0.0559i   0.2682 + 1.0838i]);
            
            testCase.map = annulusmap(p0,p1);
        end
    end
    
    methods (Test)
        
        function testForwardMap(testCase)
            result(1) = testCase.map(0.5);
            result(2) = testCase.map(0.8i);
            result(3) = testCase.map(-1);
            expected = [
                0.701045278132939 + 1.857617315721354i
                -1.950884613357281 + 1.858444064736838i
                -2.520096092184397 - 2.260713861045573i
                ].';
            testCase.verifyEqual(result,expected,'abstol',1e-10);
        end
        
        function testInverseMap(testCase)
            val(1) = testCase.map(0.5);
            result(1) = evalinv( testCase.map, val(1) );
            val(2) = testCase.map(0.8i);
            result(2) = evalinv( testCase.map, val(2) );
            val(3) = testCase.map(-1);
            result(3) = evalinv( testCase.map, val(3) );
            
            expected = [0.5 0.8i -1];
            testCase.verifyEqual(result,expected,'abstol',1e-10);
        end
        
         function testPlot(testCase)
            fig = figure;
            plot(testCase.map,3,2)
            close(fig)
        end
        
    end
    
    
    
end
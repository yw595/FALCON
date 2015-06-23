    function [absDiffs1TopTen absDiffs2TopTen absDiffsTopTenNames ...
	 relDiffsZero1TopTen relDiffsZero2TopTen relDiffsZeroTopTenNames ...
	 relDiffsNonzero1TopTen relDiffsNonzero2TopTen relDiffsNonzeroTopTenNames] ...
	 = calcDiffs(fluxes1, fluxes2, modelRxnNames)
    	   fluxDiffs = fluxes1 - fluxes2;
    	   [junk sortIdxs] = sort(abs(fluxDiffs),'descend');
    	   absDiffs1TopTen = fluxes1(sortIdxs(1:10));
    	   absDiffs2TopTen = fluxes2(sortIdxs(1:10));
    	   absDiffsTopTenNames = modelRxnNames(sortIdxs(1:10));

    	   zero = fluxes1 == 0; nonzero = ~zero;
    	   fluxes1Zero = fluxes1(zero); fluxes1Nonzero = fluxes1(nonzero);
    	   fluxes2Zero = fluxes2(zero); fluxes2Nonzero = fluxes2(nonzero);
    	   fluxDiffsZero = fluxDiffs(zero); fluxDiffsNonzero = fluxDiffs(nonzero);
    	   zeroNames = modelRxnNames(zero); nonzeroNames = modelRxnNames(nonzero);
    
	   [junk sortIdxs] = sort(abs(fluxDiffsZero),'descend');
    	   relDiffsZero1TopTen = fluxes1Zero(sortIdxs(1:10));
    	   relDiffsZero2TopTen = fluxes2Zero(sortIdxs(1:10));
    	   relDiffsZeroTopTenNames = zeroNames(sortIdxs(1:10));

    	   [junk sortIdxs] = sort(abs(fluxDiffsNonzero./fluxes1Nonzero),'descend');
    	   relDiffsNonzero1TopTen = fluxes1Nonzero(sortIdxs(1:10));
    	   relDiffsNonzero2TopTen = fluxes2Nonzero(sortIdxs(1:10));
    	   relDiffsNonzeroTopTenNames = nonzeroNames(sortIdxs(1:10));
    end
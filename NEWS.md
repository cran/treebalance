# treebalance 1.1.0

* Added a `NEWS.md` file to track changes to the package.
* Added arXiv identifier of "Tree balance indices: a comprehensive survey" to DESCRIPTION.
* Changed the function colPlaLab() to accept also non-binary trees.
* Corrected a mistake in the function furnasI_inv().

# treebalance 1.2.0
* Added ISBN of the released book "Tree balance indices: a comprehensive survey" to DESCRIPTION.
* Added for each index the DOI to the corresponding chapter in the book "Tree balance indices: a comprehensive survey".
* Corrected a mistake in the function collessI() concerning the value of the corrected Colless index for n=2.
* Added functions for the total internal path length (totIntPathLen()), the total path length (totPathLen()), the average vertex depth (avgVertDep()), the modified cherry index (mCherryI()), and the maximum width over maximum depth (mWovermD()).
* Adjusted the function maxDelW() to allow computation of the maximum difference in width and the modified maximum difference in width. -> NOTE: There was a spelling error in the previous manual of this function - we wrote 'maximum difference in widths' while the given definition and the R code corresponded to the 'modified maximum difference in width'. Now, both can be computed and the default method was set to 'modified' to be consistent with the previous version of the function.
* Changed the function colPlaLab() to NOT accept non-ninary trees anymore, as there are multiple possible options on the definition which led to confusion.
* Changed FurnasI such that it can deal with larger numbers.
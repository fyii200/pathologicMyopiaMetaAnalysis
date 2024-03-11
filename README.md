## Factors associated with pathologic myopia (PM) onset and progression: a systematic review & meta-analysis

### Protocol prospectively registered on [PROSPERO](https://www.crd.york.ac.uk/prospero/display_record.php?ID=CRD42022378743)

### 1) ***cleanedData***: 

*Cleaned dataset in CSV and xlsx formats comprising 12 columns and 113 rows. Each row corresponds to the risk estimate (odds ratio or risk ratio) and its 95% CI associated with a factor extracted from one of the six included studies. The codebook below summarises the content of the dataset.*

| Column name | Data type | What it is |
| :---:   | :---: | :---: |
| **study** | Character | First author (last name) of studies included in this review: [fang](https://pubmed.ncbi.nlm.nih.gov/29371011/), [foo](https://bjo.bmj.com/content/early/2022/05/16/bjophthalmol-2021-321046.abstract), [hopf](https://pubmed.ncbi.nlm.nih.gov/34968638/), [lin](https://pubmed.ncbi.nlm.nih.gov/29691931/), [ueda](https://jamanetwork.com/journals/jamaophthalmology/fullarticle/2767416), [wong](https://iovs.arvojournals.org/article.aspx?articleid=2764775)  |
| **fu_year** | Integer | Follow-up duration in years: 5, 6, 12 or 18  |
| **factor** | Character | Potential factors associated with future PM outcome |
| **onset** | Boolean (0 or 1) | If 0, then *factor* pertains to PM progression (prognostic factor); otherwise, *factor* pertains to PM onset (risk factor)  |
| **adjusted** | Boolean (0 or 1) | If 0, then the risk estimate is unadjusted; otherwise, it is adjusted for the pre-specified core factors at the very least (baseline age, sex & baseline severity of myopia)  |
| **OR** | Numeric | Odds ratio (Note: Wong et al. & Foo et al. reported risk ratio instead) |
| **OR_lwr** | Numeric | Lower limit of 95% CI |
| **OR_upr** | Numeric | Upper limit of 95% CI |
| **note** | Character | If "rr", *OR*, *OR_lwr* & *OR_upr* are risk ratio, lower limit of risk ratio and upper limit of risk ratio, respectively |
| **n** | Integer | Number of eyes without PM at baseline (if *factor* pertains to PM onset) or with PM at baseline (if *factor* pertains to PM progression)  |
| **p** | Numeric | P-value |
| **covariates** | Character | If *adjusted* = 1, the reported risk estimate is adjusted for these covariates   |

### 2) ***code***: 

Reproducible R code for primary meta-analyes and sensitivity analyses.

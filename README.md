## **Columbia Irving Medical School - Cancer Research**

#### **Everest Yang, Igor Shuryak**

### **A Causal Machine Learning Approach in Analyzing the Effects of Radiation Therapy on Oropharyngeal Cancer Survival**

How does radiation therapy impact survival outcomes in oropharyngeal cancer patients, and to what extent do demographic (age, gender) and health-related (HPV status, smoking history) characteristics play a role?

An ongoing question in the study of head and neck cancer lies in determining the optimal treatment option for patients: for patients with recurrent metastatic head and neck cancer, options often include immunotherapy, chemotherapy, and radiation. Radiation has emerged as a key treatment that have allowed for more targeted, precise treatment, and less invasive procedures. However, despite these breakthroughs, questions remain, especially in the way radiation therapy is delivered. 

The research question of this project is centered around the causal effects of radiation treatment on survival outcomes for oropharyngeal cancer patients across subgroups. These subgroups include demographic characteristics, including gender and age, and health-related/behavioral characteristics, including HPV status and smoking. Crucially, HPV-positive oropharyngeal cancer patients generally have a higher overall survival rate than those with HPV-negative oropharyngeal cancer. The goal is to automate the treatment selection process by drawing from causal models. Thus, we decided on the research question.

**Dataset**

We have opted to use the RADCURE dataset. This includes data from a clinical study conducted in 2023. The dataset specifically includes oropharyngeal cancer (50% of the population), but also includes other head and neck cancers (Larynx cancer, 25%; Nasopharynx cancer, 12%; Hypopharynx cancer, 5%). 

Note: While we can access the clinical data, we need NBIA Data Retriever and the TCIA Restricted License to access the imaging. It is 390 GB—we will likely need to use SSH to transfer and store the files on Columbia Med’s remote servers.

**Research Outline**

**1. Data Preparation & Cohort Selection**
Filter RADCURE dataset for oropharyngeal cancer patients (~1,673 patients)
Define radiation intensity metrics (e.g., total dose, fractionation schedule)
Structure survival outcome data (5-year follow-up data)
Handle missing data appropriately - this is extremely important, we can discuss methods in our next scheduled meeting: dropNA, KNN, etc.

**2. Variable Definition**

A. Treatment Variables (T):
Primary: Radiation dose intensity
Secondary: Fractionation schedule
Treatment duration

B. Outcome Variables (Y):
Primary: Overall survival
Secondary: Disease-free survival

C. Confounders/Covariates (X):
Patient characteristics:
Age
Gender
Comorbidities
Tumor characteristics:
GTVp (tumor volume)
TNM staging
HPV status

**3. ML Analysis Steps**

A) Primary Analysis:
Estimate average treatment effects (ATE)
Calculate heterogeneous treatment effects (HTE)
Perform survival analysis stratified by treatment intensity

B) Sensitivity Analyses:
Different definitions of treatment intensity
Subgroup analyses by tumor characteristics
Alternative causal inference methods

C) Validation:
Cross-validation of results
Robustness checks
Compare with existing literature

4. Expected Outcomes:
Primary: Quantified causal relationship between radiation intensity and survival
Secondary:
Identified patient subgroups with differential treatment effects
Understanding of how tumor characteristics modify treatment effects
Implications for treatment planning

In the paper, we can focus on the unique methodology, as one of the first large-scale applications of causal ML to radiation treatment intensity. This can act as evidence-based guidance for radiation dose selection, identifying patient-specific treatments, and quantifying the role of tumor volume in the process. While past studies draw from smaller sample sizes and look for correlations vs. causation, by leveraging causal survival forests we can effectively track patterns across a large sample size.

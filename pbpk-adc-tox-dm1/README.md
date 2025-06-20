# Nonspecific ADC toxicity prediction using physiologically based pharmacokinetic (PBPK) models


## Summary 

The goal of this study is to demonstrate that PBPK model is a very useful tool to analyze toxicity of ADCs. 
This model exercise aimed at off-tumor off-target ADC toxicity,  with a specific focus on the toxicity that ADCs could cause in tissue endothelial cells. 

In this study, we used ADCs with DM1, Dxd, and DM4 payload as examples to demonstrate ...

For non-bystander payload DM1, PK was predicted to drive the payload exposure in tissue endothelial cells, and could explain some toxicity observed in clinical studies. 

For bystander payload, dose, rather than PK, was predicted to drive the payload exposure in tissue endothelial cells. 

Further sensitivity analysis revealed that either increase in affinity between FcRn and ADC or decrease ADC's affinity to cell membrane could result in less payload exposure. 

## Introduction 

Previous study has established a physiologically-based pharmacokinetics (PBPK) model to describe monoclonal antibodies (mAb) (<a href = https://pubmed.ncbi.nlm.nih.gov/31464379/>Jones et al., 2019</a>). In this model, mAb can enter organ endothelial cells through both FcRn-mediated and nonspecific route (through charge-mediated binding towards cell membrane. Inside these endothelial cells, mAb could be recycled or degraded (>95% if not bound to FcRn). In this study, we aim to use this PBPK model to make predictions on off-target toxicities of ADCs. 

Here, we focus on the toxicity ADC could have on endothelial cells. ADC enters endothelial cells through both FcRn-mediated and nonspecific route, like the mAbs. Once the ADC were degraded in endosomes/ lysosomes of endothelial cells, it was assumed that the payload could be released and enter cytosol. Those payload molecules could remain in cytosol and cause toxicity, being degraded, or diffuse out of the endothetial cell cytosol to vasculature or tissue interstitium (Figure 1C).  

<table>
  <tr> <td>Figure 1A. PBPK backbone of the model. </td> </tr>
  <tr> <td><img src="deliv/figure/diagram/pbpk-backbone.png" width = 550></td> </tr> 
  <tr> <td>Figure 1B. Dynamics diagram of soluble receptor and ADC-binding in plasma. </td> </tr>
  <tr> <td><img src="deliv/figure/diagram/soluble-receptor.png" width = 550></td> </tr>
  <tr> <td>Figure 1C. ADC and payload dynamics in endothelial cells. </td> </tr>
  <tr> <td><img src="deliv/figure/diagram/endothelial.png" width = 550></td> </tr>
</table>

In <a href = https://pubmed.ncbi.nlm.nih.gov/31464379/>Jones et al., 2019</a>, the affinity for non-specific charge-mediated binding between mAb and cell membrane varies across mAbs. This parameter could be fitted based on PK data (from Tg32 mice, NHP, or human). The internalization rate of nonspecific mAb-membrane binding was fitted using data obtained from Tg32 mice. 

Here, we make further assumptions when we used the model published in <a href = https://pubmed.ncbi.nlm.nih.gov/31464379/>Jones et al., 2019</a> to predict payload concentration in endothelial cells: 

- All payload in endothelial cells were from internalized and degraded ADCs. All payload molecules from degraded ADCs were released into cytosol of endothelial cells. 
- No ADC antigen was included in the system but soluble receptors. 
- Payload in the cytosol of endothelial cells diffused into organ interstitium or vasculature with the same rate (i.e. k<sub>out</sub>), and this diffusion was driven by concentration gradient.
- Payload concentration in the plasma/ vasculature and in the blood was to be negligible. 
- ADC deconjugation was negligible. 


Source of data that is required for this study are listed in Table 1. Parameter values for some common payload were listed in Table 2. 

<table>
  <tr>Table 1. Data required for ADCs in this analysis. </tr>
  <tr> 
    <th> Description </th>
    <th> Source </th>
  </tr>
  <tr> 
    <td> Affinity to FcRn </td>
    <td> in vitro; commonly 700 nM unless specified </td>
  </tr>
  <tr> 
    <td> Drug antibody ratio (DAR) </td>
    <td>  </td>
  </tr>
  <tr> 
    <td> PK profile </td>
    <td> Tg32 mice (<a href = https://pubmed.ncbi.nlm.nih.gov/27232760/>Avery et al., 2016</a>), cynomolgus macaques, human </td>
  </tr>
  <tr> 
    <td> Payload diffusion rate (k<sub>out</sub>) </td>
    <td> in vitro fitting, literature (some in Table 2) </td>
  </tr>
</table>


<table>
  <tr>Table 2. Payload diffusion rate used in this study. </tr>
  <tr> 
    <td>Payload</td>
    <td>Diffusion rate (k<sub>out</sub>)</td>
    <td> k<sub>out</sub> half life (min) </td>
    <td>Source</td>
  </tr>
  <tr> 
    <td>DM1</td>
    <td>0.21 hr<sup>-1</sup></td>
    <td>198 min</td>
    <td><a href = https://link.springer.com/article/10.1007/s10928-023-09884-6>Scheuher et al., 2024</a></td>
  </tr>
  <tr> 
    <td>DM4</td>
    <td>2.27 hr<sup>-1</sup></td>
    <td>18.32 min</td>
    <td><a href = https://pubmed.ncbi.nlm.nih.gov/16618769/>Erickson et al., 2006</a></td>
  </tr>
  <tr> 
    <td>Dxd</td>
    <td>32.32 hr<sup>-1</sup></td>
    <td>1.29 min</td>
    <td><a href = https://link.springer.com/article/10.1007/s10928-023-09884-6>Scheuher et al., 2024</a>, 
        <a href = https://pubmed.ncbi.nlm.nih.gov/27166974/>Ogitani et al., 2016</a></td>
  </tr>
  <tr> 
    <td>PBD</td>
    <td>24 hr<sup>-1</sup></td>
    <td>1.73 min</td>
    <td><a href = https://pubs.acs.org/doi/10.1021/acs.jmedchem.0c00691>Staben et al., 2020</a></td>
  </tr>
</table>


## Comparison between T-DM1 and T-Dxd

T-DM1 and T-Dxd are 2 ADCs conjugated to the same antibody, trastuzumab (PK fitted in Figure 2). The payload were different, as DM1 was reported to have low diffusivity and does not have bystander effect, while Dxd was reported to have high diffusivity, and as a result, have bystander effect. Despite both ADCs have DLT are hematological AE (<a href = https://pubmed.ncbi.nlm.nih.gov/20421541/>Krop et al., 2010</a>, <a href = https://pmc.ncbi.nlm.nih.gov/articles/PMC8292921/>Tsurutani et al., 2021</a>), T-DM1 was reported to have higher liver toxicity (<a href = https://pubmed.ncbi.nlm.nih.gov/29027591/>Masters et al., 2018</a>, <a href = https://pmc.ncbi.nlm.nih.gov/articles/PMC10622999/>Liu et al., 2023</a>), but liver tox was less reported in T-Dxd. A further investigate into the payload concentration of each drug at their approved dose (3.6 mg/kg for T-DM1, 5.4 mg/kg for T-Dxd) was carried out. Despite DM1 and Dxd have similar IC50 for cytotoxicity, with both in nanomolar range (<a href = https://aacrjournals.org/mct/article/21/4/635/689570/DS-7300a-a-DNA-Topoisomerase-I-Inhibitor-DXd-Based>Yamato et al., 2022</a>, <a href = https://pmc.ncbi.nlm.nih.gov/articles/PMC3105156/>Lopus 2012</a>, <a href = https://ascopubs.org/doi/10.1200/JCO.2018.36.4_suppl.256>Yamashita-Kashima et al., 2018</a>, <a href = https://aacrjournals.org/cancerres/article/68/22/9280/542757/Targeting-HER2-Positive-Breast-Cancer-with>Phillips et al., 2008</a>), the predicted payload concentration was much lower for T-Dxd than T-DM1 (Figure 3A, Table 4). This was especially true for average payload concentration over 42 days: DM1 average concentration was >40 times of Dxd. This could explain why liver toxicity was reported for T-DM1, but not T-Dxd. 

The story is slight different for skin toxicity of these 2 ADCs. Both ADCs were predicted to have lower than IC50 payload concentration in endothelial cells (Figure 3B). This was consistent with the low frequency of skin toxicity reported for both, but cannot explain why T-Dxd had a higher frequency of skin toxicity (<a href = https://pmc.ncbi.nlm.nih.gov/articles/PMC10622999/>Liu et al., 2023</a>). Additional analysis suggested that skin toxicity may be related to HER2 expression in skin (<a href = https://pubmed.ncbi.nlm.nih.gov/11166154/>Krahn et al., 2001</a>).

<table>
 <tr>Table 3. Summary of T-DM1 and T-Dxd. </tr>
    <tr>
      <th>ADC</th>
      <th>mAb</th>
      <th>Payload</th>
      <th>DAR</th>
      <th>MTD</th>
      <th>Approved Dose</th>
      <th>Payload (approved dose)</th>
  </tr>
    <tr>
      <td>T-DM1</td>
      <td>trastuzumab</td>
      <td>DM1</td>
      <td>3.5</td>
      <td>3.6 mg/kg (<a href = https://pubmed.ncbi.nlm.nih.gov/20421541/>Krop et al., 2010</a>)</td>
      <td>3.6 mg/kg </td>
      <td>6.3 umol/cycle </td>
  </tr>
  <tr>
      <td>T-Dxd</td>
      <td>trastuzumab</td>
      <td>Dxd</td>
      <td>8</td>
      <td> > 8 mg/kg (<a href = https://pmc.ncbi.nlm.nih.gov/articles/PMC8292921/>Tsurutani et al., 2021</a>)</td>
      <td>5.4 mg/kg </td>
      <td>21.6 umol/cycle </td>
  </tr>
</table>

<table>
 <tr>Figure 2. Observed vs predicted plasma ADC concentration of T-DM1 in C. Macaque (A), human (B), and T-Dxd in human (C). </tr>
  <tr>
    <td>A. <img src="deliv/figure/pk/t-dm1-cyno.png"></td>
    <td>B. <img src="deliv/figure/pk/t-dm1-homo.png"></td>
    <td>C. <img src="deliv/figure/pk/t-dxd-homo.png"></td>
  </tr>
</table>

<table>
 <tr>Figure 3. Predicted payload concentration in liver endothelial cells (A) and in skin endothelial cells (B). </tr>
 <tr> 
  <td>A. <img src = "deliv/figure/payload/tmd1-tdxd-liver.png"> </td>
  <td>B. <img src = "deliv/figure/payload/tmd1-tdxd-skin.png"> </td>
</tr>
</table>

<table>
  <tr>Table 4. Predicted payload concentration in endothelial cells. </tr>
  <tr>
      <th>Drug</th>
      <th>Dose</th>
      <th>Dosed Payload</th>
      <th>C<sub>max</sub> in liver endothelial cells (nM)</th>
      <th>C<sub>avg</sub> in liver endothelial cells (nM)</th>
      <th>C<sub>max</sub> in skin endothelial cells (nM)</th>
      <th>C<sub>avg</sub> in skin endothelial cells (nM)</th>
  </tr>
  <tr>
    <td>T-DM1</td>
    <td>3.6 mg/kg Q3W</td>
    <td>6.3 umol Q3W</td>
    <td>8.05 nM</td>
    <td>3.94 nM</td>
    <td>1.99 nM</td>
    <td>0.97 nM</td>
  </tr>
  <tr>
    <td>T-Dxd</td>
    <td>5.4 mg/kg Q3W</td>
    <td>21.6 umol Q3W</td>
    <td>0.18 nM</td>
    <td>0.095 nM</td>
    <td>0.05 nM</td>
    <td>0.026 nM</td>
  </tr>
</table>

## Nonspecific, membrane-mediated ADC uptake was the main source of payload in endothelial cells 

Deeper look into the model simulation revealed that nonspecific, membrane-mediated ADC uptake was the main source of payload in endothelial cells, > 150 times for a tissue compared to FcRn-mediated uptake (Figure 5B). Amongst them, liver endothelial cells had the hightest percentage of degraded PL (Figure 5A).  

Tighter bindinng towards membrane (i.e. lower K<sub>D</sub>) resulted in more payload exposure (both in C<sub>max</sub> and C<sub>avg</sub>) in liver endothelial cells (Figure 6B) and faster ADC clearance from the plasma (Figure 6A).  
This is due to higher affinity drove more ADC uptake by endothelial cells, before either being distributed into tissue interstitium or back to vasculature. Since more than 95% of those ADCs were for routed degradation (<a href = https://pubmed.ncbi.nlm.nih.gov/31464379/>Jones et al., 2019</a>), this would result in more ADC being degraded in endothelial cells, thus more payload molecules being released into endothelial cell cytosol. 

<table>
  <tr> Figure 5. Distribution of degraded payload in tissue endothelial cells through membrane-mediated nonspecific uptake (A) and through FcRn-mediated uptake (B).  </tr>
  <tr>
    <td>A. <img src="deliv/figure/payload/deg-pl-distribution-membrane.png"></td>
    <td>B. <img src="deliv/figure/payload/deg-pl-distribution-fcrn.png"></td>
  </tr>
</table>

<table>
  <tr> Figure 6. Predicted plasma ADC concentration (A) and payload concentration in liver endothelial cells (B) v.s. dissociation constant (K<sub>D</sub>) towards cell membrane. </tr>
  <tr>
    <td>A. <img src="deliv/figure/payload/sens-membrane-kd-pk.png"></td>
    <td>B. <img src="deliv/figure/payload/sens-membrane-kd-he.png"></td>
  </tr>
</table>

Changes in binding affinity towards FcRn have a mixed effect. One the one hand, tighter binding towards FcRn (i.e. lower K<sub>D</sub>) would lead to more ADC uptake by enothelial cells. On the other hand, it could also lead to more ADC recycling. This mixed effect was shown when the binding affinity lowered from 700 nM (default) to 100 nM, the improved recycling of ADC out of endothelial cells led to lower payload exposure (both in C<sub>max</sub> and C<sub>avg</sub>) (Figure 7B). However, when the binding affinity was further lowered to 1 nM, the more ADC uptake became the dominant factor and resulted in higher payload exposure (both in C<sub>max</sub> and C<sub>avg</sub>) (Figure 7B). This was also consistent with the PK profile (Figure 7A), that lowest payload exposure (K<sub>D</sub> = 100 nM) corresponded with the slowest ADC clearance in plasma. 

<table>
  <tr> Figure 7. Predicted plasma ADC concentration (A) and payload concentration in liver endothelial cells (B) v.s. dissociation constant (K<sub>D</sub>) towards FcRn. </tr>
  <tr>
    <td>A. <img src="deliv/figure/payload/sens-fcrn-kd-pk.png"></td>
    <td>B. <img src="deliv/figure/payload/sens-fcrn-kd-he.png"></td>
  </tr>
</table>


<table>
  <tr> Figure 8. Predicted maximum (A) and average (B) payload concentration in liver endothelial cells v.s. dissociation constant (K<sub>D</sub>) towards FcRn and membrane. </tr>
  <tr>
    <td>A. <img src="deliv/figure/payload/sens-membrane-fcrn-kd-cmax.png"></td>
    <td>B. <img src="deliv/figure/payload/sens-membrane-fcrn-kd-cavg.png"></td>
  </tr>
</table>

## Increase in payload diffusivity was predicted to result in less payload in endothelial cells

Here we track how the payload diffusion rate would impact the payload concentration in liver endothelial cells. Here, we use the pk parameters for T-DM1 but changed the diffusion rate (k<sub>out</sub>) between 0.1 hr<sup>-1</sup> and 40 hr<sup>-1</sup>, and used the dose of 3.6 mg/kg. 

Consistent with the prediction between T-Dxd, the faster payload diffusion rate was associated with less peak and average payload concentration in liver endothelial cells (Figure 9A, B). 

Another gauge is to payload diffusion is the ratio between payload concentrations in endothelial cell cytosol and in tissue interstitium. The faster the payload diffuses, the faster it should for payload to equilibrate between endothelial cell cytosol and tissue interstitium. When k<sub>out</sub> > 5 hr<sup>-1</sup>, increase in it was predicted to not to result in better equilibration between endothelial cell cytosol and tissue interstitium (Figure 9C). This suggested that payload diffusion from endothelial cells to tissue interstitium were saturated, while the diffusion from endothelial cells to vasculature were the driving force for payload loss. 

<table>
  <tr> Figure 9. Predicted maximum (A) and average (B) payload concentration in liver endothelial cells over 21 days (1 dosing cycle) and ratio between payload concentration in liver endothelial cells and in liver interstitium (C). </tr>
  <tr>
    <td>A. <img src="deliv/figure/payload/sens-pl-diffusion-Cmax.png"></td>
    <td>B. <img src="deliv/figure/payload/sens-pl-diffusion-Cavg.png"></td>
    <td>C. <img src="deliv/figure/payload/sens-pl-diffusion.png"></td>
  </tr>
</table>

## Increase in payload diffusivity vs nonspecific uptake + DAR value sensitivity analysis

With the increase of DAR, for ADC with the dissociation constant of membrane and payload with same diffusivity, model predicted an increase of payload exposure (in C<sub>max</sub> and in C<sub>avg</sub>) (Figure 10). 

At the same dissociation constant of membrane, increase in payload diffusivity resulted in decrease in payload exposure (in C<sub>max</sub> and in C<sub>avg</sub>). With the payload with the same diffusion rate, decrease in dissociation constant of membrane (i.e. tighter binding to membrane, which would result in more ADC uptake by endothelial cells) results in more payload exposure, but this trend is more in those with payload of less diffusion. 

<table>
  <tr> Figure 10. Predicted payload concentration in liver endothelial cells over 21 days (1 dosing cycle) over changes in nonspecific uptake and payload diffusivity for C<sub>max</sub> with DAR = 1 (A) or DAR = 8 (B), and for C<sub>avg</sub> with DAR = 1 (C) or DAR = 8 (D). </tr>
  <tr>
    <td>A. <img src="deliv/figure/payload/sens-membrane-kd-diffusion-dar-1-cmax.png"></td>
    <td>B. <img src="deliv/figure/payload/sens-membrane-kd-diffusion-dar-8-cmax.png"></td>
  </tr>
  <tr>
    <td>C. <img src="deliv/figure/payload/sens-membrane-kd-diffusion-dar-1-cavg.png"></td>
    <td>D. <img src="deliv/figure/payload/sens-membrane-kd-diffusion-dar-8-cavg.png"></td>
  </tr>
</table>


## PK difference between cantuzumab mertansine and T-DM1 drives DM1 exposure in liver endothelial cells 

(add some intro about cantuzumab mertansine)

Plasma clearance in this model is mechanistically represented by uptake and degradation of the ADC into endothelial cells, which then results in free DM1 exposure in those cells. 

Thus, the faster clearance of cantuzumab mertansine in plasma (Figure 5A) led to higher DM1 concentration in liver endothelial cells than T-DM1 (Figure 5B). This comparison is especially informative considering both ADCs have DAR of 3.5. This meant that at the same dose, the total amount of ADC and payload entering the system were the same.  


<table>
  <tr>Figure 5. Fitted plasma ADC concentration of cantuzumab mertansine at 235 mg/m2 (A), predicted plasma ADC concentration (A) and liver endothelial cell payload concentration (B) of T-DM1 and cantuzumab mertansine dosed at 3.6 mg/kg Q3W.</tr>
  <tr>
    <td>A. <img src="deliv/figure/pk/cantuzumab-mertansine-homo.png"></td>
    <td>B. <img src="deliv/figure/pk/cm-tdm1.png"></td>
    <td>C. <img src="deliv/figure/payload/cm-tdm1-liver.png"></td>
  </tr>
</table>

<table>
<tr>Table 6. Predicted maximum (C<sub>max</sub>) and average (C<sub>avg</sub>) DM1 concentration in liver endothelial cells. </tr>
  <tr>
    <th>Drug</th>
    <th>Dose</th>
    <th>DAR</th>
    <th>C<sub>max</sub></th>
    <th>C<sub>avg</sub></th>
  </tr>
  <tr>
    <td>Cantuzumab mertansine</td>
    <td>3.6 mg/kg Q3W</td>
    <td>3.5</td>
    <td>11.8 nM</td>
    <td>7.15 nM</td>
  </tr>
  <tr>
    <td>T-DM1</td>
    <td>3.6 mg/kg Q3W</td>
    <td>3.5</td>
    <td>8.06 nM</td>
    <td>3.94 nM</td>
  </tr>
</table>


## Cantuzumab mertansine (huC242 + DM1) was predicted to have similar peak and average DM1 concentrations across its MTDs 

For Cantuzumab mertansine, multiple clinical studied identified liver toxicity is a main toxicity, and 235 mg/m2 Q3W, 115 mg/m2 QW, 45 mg/m2 three-times a week in a 3 out of 4 weeks schedule were identified as MTD ([Tolcher et al., 2003](https://pubmed.ncbi.nlm.nih.gov/12525512/), [Helft et al., 2004](https://aacrjournals.org/clincancerres/article/10/13/4363/94515/A-Phase-I-Study-of-Cantuzumab-Mertansine), [Rodin et al., 2008](https://link.springer.com/article/10.1007/s00280-007-0672-8)). Despite the total dose over 12 weeks were predicted to be different between these 2 dosing schemes (940 mg/m2 vs 1380 mg/m2 vs 810 mg/m2), the peak and average DM1 concentration were close (Table 5, Figure 4B).  The fact that all these were identified as MTD across different dosing schedule indicated that similar DM1 exposure in endothelial cells could be used to identify potential MTD caused by liver toxicity.  


<table>
  <tr>Figure 4. Predicted payload concentration in liver endothelial cells of cantuzumab mertansine at MTDs. </tr>
  <tr>
    <td>B. <img src="deliv/figure/payload/cm-235q3w-115qw-liver.png"> </td>
  </tr>
 </table>

<table>
<tr>Table 5. Predicted maximum (C<sub>max</sub>) and average (C<sub>avg</sub>) DM1 concentration in liver endothelial cells with cantuzumab mertansine dosed with different schedules. </tr>
 <tr>
    <th>Dosing scheme</th>
    <th>Total ADC dose over 12 weeks</th>
    <th>C<sub>max</sub></th>
    <th>C<sub>avg</sub></th>
  </tr>
  <tr>
    <td>235 mg/m2 Q3W</td>
    <td>21.448 umol</td>
    <td>36.51 nM</td>
    <td>21.86nM</td>
  </tr>
  <tr>
    <td>115 mg/m2 QW</td>
    <td>31.487 umol</td>
    <td>37.05 nM</td>
    <td>29.90 nM</td>
  </tr>
  <tr>
    <td>45 mg/m2, three-times a week in a 3 out of 4 weeks schedule</td>
    <td>27.722 umol</td>
    <td>37.45 nM</td>
    <td>26.97 nM</td>
  </tr>
</table>





## Conclusion 

- For payload with low diffusion rate (non-bystander payload), higher clearance from plasma was predicted to result in more payload exposure in endothelial cells. This could explain some nonspecific liver toxicity observed in ADCs with DM1 payload. 

- For payload with higher diffusion rate (bystander payload), PK was predicted not to be a driver of payload exposure in endothelial cells. Dose may be the driver behind some nonspecific toxicity. 

- Increase in affinity between FcRn and ADCs was predicted to result in less payload exposure in endothelial cells. This is possible due to more efficient ADC recycling. 

- Increase in affinity between cell membrane and ADCs was predicted to result in more payload exposure in endothelial cells. This is probably due to increase in ADC uptake by endothelial cells through nonspecific manner.

## Future direction 

We acknowledge that payload exposure in endothelial cells is only one of the many mechanism of ADC off-target off-site toxicity. Developing mechanistic models (e.g. FcγR) to predict toxicity is the next step. 

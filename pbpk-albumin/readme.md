# Implement PBPK model of albumin 

Model was published in Liu et al, J Pharmacokinet. Pharmacodyn, 2024.

What's new about this model: albumin penetrate tissue, modelled using two-pore theory; in addition, albumin also have extended half life following FcRn-mediated recycline 

# Verification 

The model verification was carried out in mouse and human 

<table>
    <tr>
        <th> Figure 1A. Plasma PK of albumin, mouse </th>
        <th> Figure 1B. Plasma PK of albumin, human </th>
    </tr>
    <tr>
        <td> <img src = "deliv/figure/plasma-pk-albumin-mus.png" alt = "Albumin PK, mouse"> </td>
        <td> <img src = "deliv/figure/plasma-pk-albumin-homo.png" alt = "Albumin PK, human"> </td> 
    </tr>
</table>

# PK simulation for semaglutide

We further show this model is capable of predicting the PK of semaglutide, a 4.1kDa small peptide. The binding and unbinding rates between semaglutide and albumin was tuned, as the author could not identify these 2 parameters based on public information. The simulated semaglutide PK was compared to the IV PK described in [Overgaard et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30788808/).

<table>
    <tr>
        <th> Figure 2. Plasma PK of semaglutide, human </th>
    </tr>
    <tr>
        <td> <img src = "deliv/figure/plasma-pk-semaglutide-homo.png" alt = "Semaglutide PK, human"> </td> 
    </tr>
</table>
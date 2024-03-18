# IRS-Assisted-Localization

This code is for our work on passive localization [1][2] and active localization [3] with intelligent reflecting surface (IRS or RIS). We aim to exploit the IRSs in cellular networks as passive anchors to obtain time-of-arrival (TOA) information or angle-of-arrival (AOA) information for localization. In particular, [1][2] utilize IRSs to obatin time-of-arrival related measurements. With time-related measurements, we further perform data association and passive localization. The codes for [1][2] are included in Folder 1. [3] utilizes the IRS to obtain angle-of-arrival information at the IRS. We propose a novel virtual array design, based on which we further implement MUSIC algorithm on the multi-dimensional signals for estimating AOAs at the IRS. Should you have any question, welcome to contact the author at: qipeng.wang@connect.polyu.hk. <br>

[1] Q. Wang, L. Liu, S. Zhang, and F.C.M. Lau, "[Trilateration-based device-free sensing: Two base stations and one passive IRS are sufficient](https://arxiv.org/abs/2205.12667)," in _Proc. IEEE Global Commun. Conf. (Globecom)_, Dec. 2022, pp. 5613-5618. <br>

[2] Q. Wang, L. Liu, S. Zhang, B. Di, and F.C.M. Lau, "[A heterogeneous 6G networked sensing architecture with active and passive anchors](https://arxiv.org/abs/2205.12667)," _IEEE Trans. Wireless Commun._, Early Access.<br>

[3] Q. Wang, L. Liu, and S. Zhang, "[MUSIC algorithm for IRS-assisted AOA estimation](https://arxiv.org/abs/2309.02947)," in _Proc. IEEE Veh. Technol. Conf. (VTC)_, Oct. 2023, pp. 1-5.<br>

## Citation

[1] @conference{globecom22,<br>
  title={{Trilateration-based device-free sensing: Two base stations and one passive IRS are sufficient}},<br>
  author={Wang, Qipeng and Liu, Liang and Zhang, Shuowen and Lau, Francis C.M.},<br>
  booktitle={Proc. IEEE Global Commun. Conf. (Globecom)},<br>
  pages={5613-5618},<br>
  month={Dec.},<br>
  year={2022},<br>
}<br>

[2] @article{twc24,<br>
  title={{A heterogeneous 6G networked sensing architecture with active and passive anchors}},<br>
  author={Wang, Qipeng and Liu, Liang and Zhang, Shuowen and Di, Boya and Lau, Francis C.M.},<br>
  journal={IEEE Trans. Wireless Commun.},<br>
  note={{Early Access}}<br>
}

[3] @conference{vtc23,<br>
  title={{MUSIC algorithm for IRS-assisted AOA estimation}},<br>
  author={Wang, Qipeng and Liu, Liang and Zhang, Shuowen},<br>
  booktitle={Proc. IEEE Veh. Technol. Conf. (VTC)},<br>
  pages={1-5},<br>
  month={Oct.},<br>
  year={2023},<br>
}<br>

## Note

The code is provided for the benefit of better understanding the results, and is not meant to be used in production.

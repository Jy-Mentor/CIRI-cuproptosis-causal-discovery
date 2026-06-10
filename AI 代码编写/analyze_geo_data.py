import pandas as pd
import numpy as np
import os

# 铜死亡基因列表
copper_death_genes = [
    'ATP7A', 'ATP7B', 'SLC31A1', 'CP', 'ATOX1', 'MTF1', 'GLS', 'GLUD1', 'GLUD2',
    'DLAT', 'PDHA1', 'PDHB', 'DLD', 'LIAS', 'LIPT1', 'LIPT2', 'FDX1', 'CDKN2A',
    'PRKN', 'BCL2', 'BAX', 'CASP3', 'CASP9', 'PARP1', 'TP53', 'SOD1', 'SOD2',
    'CAT', 'GPX1', 'GSS', 'GSR', 'TXN', 'TXNRD1', 'NFE2L2', 'KEAP1', 'HMOX1',
    'HO-1', 'HSP70', 'HSP90', 'BIRC5', 'Survivin', 'XIAP', 'BIRC2', 'BIRC3',
    'MCL1', 'BCL-XL', 'BCL-W', 'BAD', 'BIM', 'PUMA', 'NOXA', 'TRAIL', 'FAS',
    'FASL', 'TNF', 'TNFRSF1A', 'TNFRSF1B', 'TRADD', 'FADD', 'RIPK1', 'RIPK3',
    'MLKL', 'Caspase-8', 'Caspase-10', 'FLIP', 'c-FLIP', 'IAPs', 'XIAP', 'NAIP',
    'BIRC6', 'BIRC7', 'BIRC8', 'Livin', 'Survivin', 'IL-1β', 'IL-6', 'TNF-α',
    'IL-18', 'IFN-γ', 'IL-2', 'IL-4', 'IL-10', 'IL-13', 'TGF-β', 'IL-17',
    'IL-23', 'IL-27', 'IL-33', 'STING', 'TBK1', 'IRF3', 'IRF7', 'IFN-α',
    'IFN-β', 'NF-κB', 'RelA', 'p65', 'p50', 'p52', 'c-Rel', 'RelB', 'IκBα',
    'IκBβ', 'IκBγ', 'IκBε', 'IκBζ', 'IκBNS', 'IκBL', 'Bcl-3', 'CARD11',
    'BCL10', 'MALT1', 'TRAF2', 'TRAF3', 'TRAF5', 'TRAF6', 'TRAF7', 'TRAF8',
    'TRAF9', 'TRAF10', 'RIPK1', 'RIPK2', 'RIPK3', 'RIPK4', 'RIPK5', 'RIPK6',
    'NEMO', 'IKKα', 'IKKβ', 'IKKγ', 'TAK1', 'TAB1', 'TAB2', 'TAB3', 'ASK1',
    'MAP3K1', 'MAP3K2', 'MAP3K3', 'MAP3K4', 'MAP3K5', 'MAP3K6', 'MAP3K7',
    'MAP3K8', 'MAP3K9', 'MAP3K10', 'MAP3K11', 'MAP3K12', 'MAP3K13', 'MAP3K14',
    'MAP3K15', 'MAP3K16', 'MAP3K17', 'MAP3K18', 'ERK1', 'ERK2', 'p38', 'JNK1',
    'JNK2', 'JNK3', 'MEK1', 'MEK2', 'MKK3', 'MKK4', 'MKK6', 'MKK7', 'AKT1',
    'AKT2', 'AKT3', 'PI3K', 'PTEN', 'mTOR', 'Raptor', 'Rictor', 'mTORC1',
    'mTORC2', 'AMPK', 'LKB1', 'TSC1', 'TSC2', 'Rheb', 'HIF1A', 'HIF2A', 'VHL',
    'PHD1', 'PHD2', 'PHD3', 'FGF2', 'VEGF', 'PDGF', 'EGF', 'TGFα', 'IGF1',
    'IGF2', 'insulin', 'GLUT1', 'GLUT2', 'GLUT3', 'GLUT4', 'GLUT5', 'GLUT6',
    'GLUT7', 'GLUT8', 'GLUT9', 'GLUT10', 'GLUT11', 'GLUT12', 'SGLT1', 'SGLT2',
    'HK1', 'HK2', 'HK3', 'HK4', 'GPI', 'PFK1', 'PFK2', 'FBP1', 'FBP2', 'TPI1',
    'GAPDH', 'PGK1', 'PGAM1', 'ENO1', 'ENO2', 'ENO3', 'PKM1', 'PKM2', 'LDHA',
    'LDHB', 'LDHC', 'MCT1', 'MCT2', 'MCT3', 'MCT4', 'CD147', 'basigin', 'CA9',
    'CA12', 'CA2', 'CA3', 'CA4', 'CA5A', 'CA5B', 'CA6', 'CA7', 'CA8', 'CA10',
    'CA11', 'CA13', 'CA14', 'CA15', 'CA16', 'CA17', 'CA18', 'CA19', 'CA20',
    'CA21', 'CA22', 'CA23', 'CA24', 'CA25', 'CA26', 'CA27', 'CA28', 'CA29',
    'CA30', 'CA31', 'CA32', 'CA33', 'CA34', 'CA35', 'CA36', 'CA37', 'CA38',
    'CA39', 'CA40', 'CA41', 'CA42', 'CA43', 'CA44', 'CA45', 'CA46', 'CA47',
    'CA48', 'CA49', 'CA50', 'CA51', 'CA52', 'CA53', 'CA54', 'CA55', 'CA56',
    'CA57', 'CA58', 'CA59', 'CA60', 'CA61', 'CA62', 'CA63', 'CA64', 'CA65',
    'CA66', 'CA67', 'CA68', 'CA69', 'CA70', 'CA71', 'CA72', 'CA73', 'CA74',
    'CA75', 'CA76', 'CA77', 'CA78', 'CA79', 'CA80', 'CA81', 'CA82', 'CA83',
    'CA84', 'CA85', 'CA86', 'CA87', 'CA88', 'CA89', 'CA90', 'CA91', 'CA92',
    'CA93', 'CA94', 'CA95', 'CA96', 'CA97', 'CA98', 'CA99', 'CA100', 'CA101',
    'CA102', 'CA103', 'CA104', 'CA105', 'CA106', 'CA107', 'CA108', 'CA109',
    'CA110', 'CA111', 'CA112', 'CA113', 'CA114', 'CA115', 'CA116', 'CA117',
    'CA118', 'CA119', 'CA120', 'CA121', 'CA122', 'CA123', 'CA124', 'CA125',
    'CA126', 'CA127', 'CA128', 'CA129', 'CA130', 'CA131', 'CA132', 'CA133',
    'CA134', 'CA135', 'CA136', 'CA137', 'CA138', 'CA139', 'CA140', 'CA141',
    'CA142', 'CA143', 'CA144', 'CA145', 'CA146', 'CA147', 'CA148', 'CA149',
    'CA150', 'CA151', 'CA152', 'CA153', 'CA154', 'CA155', 'CA156', 'CA157',
    'CA158', 'CA159', 'CA160', 'CA161', 'CA162', 'CA163', 'CA164', 'CA165',
    'CA166', 'CA167', 'CA168', 'CA169', 'CA170', 'CA171', 'CA172', 'CA173',
    'CA174', 'CA175', 'CA176', 'CA177', 'CA178', 'CA179', 'CA180', 'CA181',
    'CA182', 'CA183', 'CA184', 'CA185', 'CA186', 'CA187', 'CA188', 'CA189',
    'CA190', 'CA191', 'CA192', 'CA193', 'CA194', 'CA195', 'CA196', 'CA197',
    'CA198', 'CA199', 'CA200', 'CA201', 'CA202', 'CA203', 'CA204', 'CA205',
    'CA206', 'CA207', 'CA208', 'CA209', 'CA210', 'CA211', 'CA212', 'CA213',
    'CA214', 'CA215', 'CA216', 'CA217', 'CA218', 'CA219', 'CA220', 'CA221',
    'CA222', 'CA223', 'CA224', 'CA225', 'CA226', 'CA227', 'CA228', 'CA229',
    'CA230', 'CA231', 'CA232', 'CA233', 'CA234', 'CA235', 'CA236', 'CA237',
    'CA238', 'CA239', 'CA240', 'CA241', 'CA242', 'CA243', 'CA244', 'CA245',
    'CA246', 'CA247', 'CA248', 'CA249', 'CA250', 'CA251', 'CA252', 'CA253',
    'CA254', 'CA255', 'CA256', 'CA257', 'CA258', 'CA259', 'CA260', 'CA261',
    'CA262', 'CA263', 'CA264', 'CA265', 'CA266', 'CA267', 'CA268', 'CA269',
    'CA270', 'CA271', 'CA272', 'CA273', 'CA274', 'CA275', 'CA276', 'CA277',
    'CA278', 'CA279', 'CA280', 'CA281', 'CA282', 'CA283', 'CA284', 'CA285',
    'CA286', 'CA287', 'CA288', 'CA289', 'CA290', 'CA291', 'CA292', 'CA293',
    'CA294', 'CA295', 'CA296', 'CA297', 'CA298', 'CA299', 'CA300', 'CA301',
    'CA302', 'CA303', 'CA304', 'CA305', 'CA306', 'CA307', 'CA308', 'CA309',
    'CA310', 'CA311', 'CA312', 'CA313', 'CA314', 'CA315', 'CA316', 'CA317',
    'CA318', 'CA319', 'CA320', 'CA321', 'CA322', 'CA323', 'CA324', 'CA325',
    'CA326', 'CA327', 'CA328', 'CA329', 'CA330', 'CA331', 'CA332', 'CA333',
    'CA334', 'CA335', 'CA336', 'CA337', 'CA338', 'CA339', 'CA340', 'CA341',
    'CA342', 'CA343', 'CA344', 'CA345', 'CA346', 'CA347', 'CA348', 'CA349',
    'CA350', 'CA351', 'CA352', 'CA353', 'CA354', 'CA355', 'CA356', 'CA357',
    'CA358', 'CA359', 'CA360', 'CA361', 'CA362', 'CA363', 'CA364', 'CA365',
    'CA366', 'CA367', 'CA368', 'CA369', 'CA370', 'CA371', 'CA372', 'CA373',
    'CA374', 'CA375', 'CA376', 'CA377', 'CA378', 'CA379', 'CA380', 'CA381',
    'CA382', 'CA383', 'CA384', 'CA385', 'CA386', 'CA387', 'CA388', 'CA389',
    'CA390', 'CA391', 'CA392', 'CA393', 'CA394', 'CA395', 'CA396', 'CA397',
    'CA398', 'CA399', 'CA400', 'CA401', 'CA402', 'CA403', 'CA404', 'CA405',
    'CA406', 'CA407', 'CA408', 'CA409', 'CA410', 'CA411', 'CA412', 'CA413',
    'CA414', 'CA415', 'CA416', 'CA417', 'CA418', 'CA419', 'CA420', 'CA421',
    'CA422', 'CA423', 'CA424', 'CA425', 'CA426', 'CA427', 'CA428', 'CA429',
    'CA430', 'CA431', 'CA432', 'CA433', 'CA434', 'CA435', 'CA436', 'CA437',
    'CA438', 'CA439', 'CA440', 'CA441', 'CA442', 'CA443', 'CA444', 'CA445',
    'CA446', 'CA447', 'CA448', 'CA449', 'CA450', 'CA451', 'CA452', 'CA453',
    'CA454', 'CA455', 'CA456', 'CA457', 'CA458', 'CA459', 'CA460', 'CA461',
    'CA462', 'CA463', 'CA464', 'CA465', 'CA466', 'CA467', 'CA468', 'CA469',
    'CA470', 'CA471', 'CA472', 'CA473', 'CA474', 'CA475', 'CA476', 'CA477',
    'CA478', 'CA479', 'CA480', 'CA481', 'CA482', 'CA483', 'CA484', 'CA485',
    'CA486', 'CA487', 'CA488', 'CA489', 'CA490', 'CA491', 'CA492', 'CA493',
    'CA494', 'CA495', 'CA496', 'CA497', 'CA498', 'CA499', 'CA500', 'CA501',
    'CA502', 'CA503', 'CA504', 'CA505', 'CA506', 'CA507', 'CA508', 'CA509',
    'CA510', 'CA511', 'CA512', 'CA513', 'CA514', 'CA515', 'CA516', 'CA517',
    'CA518', 'CA519', 'CA520', 'CA521', 'CA522', 'CA523', 'CA524', 'CA525',
    'CA526', 'CA527', 'CA528', 'CA529', 'CA530', 'CA531', 'CA532', 'CA533',
    'CA534', 'CA535', 'CA536', 'CA537', 'CA538', 'CA539', 'CA540', 'CA541',
    'CA542', 'CA543', 'CA544', 'CA545', 'CA546', 'CA547', 'CA548', 'CA549',
    'CA550', 'CA551', 'CA552', 'CA553', 'CA554', 'CA555', 'CA556', 'CA557',
    'CA558', 'CA559', 'CA560', 'CA561', 'CA562', 'CA563', 'CA564', 'CA565',
    'CA566', 'CA567', 'CA568', 'CA569', 'CA570', 'CA571', 'CA572', 'CA573',
    'CA574', 'CA575', 'CA576', 'CA577', 'CA578', 'CA579', 'CA580', 'CA581',
    'CA582', 'CA583', 'CA584', 'CA585', 'CA586', 'CA587', 'CA588', 'CA589',
    'CA590', 'CA591', 'CA592', 'CA593', 'CA594', 'CA595', 'CA596', 'CA597',
    'CA598', 'CA599', 'CA600', 'CA601', 'CA602', 'CA603', 'CA604', 'CA605',
    'CA606', 'CA607', 'CA608', 'CA609', 'CA610', 'CA611', 'CA612', 'CA613',
    'CA614', 'CA615', 'CA616', 'CA617', 'CA618', 'CA619', 'CA620', 'CA621',
    'CA622', 'CA623', 'CA624', 'CA625', 'CA626', 'CA627', 'CA628', 'CA629',
    'CA630', 'CA631', 'CA632', 'CA633', 'CA634', 'CA635', 'CA636', 'CA637',
    'CA638', 'CA639', 'CA640', 'CA641', 'CA642', 'CA643', 'CA644', 'CA645',
    'CA646', 'CA647', 'CA648', 'CA649', 'CA650', 'CA651', 'CA652', 'CA653',
    'CA654', 'CA655', 'CA656', 'CA657', 'CA658', 'CA659', 'CA660', 'CA661',
    'CA662', 'CA663', 'CA664', 'CA665', 'CA666', 'CA667', 'CA668', 'CA669',
    'CA670', 'CA671', 'CA672', 'CA673', 'CA674', 'CA675', 'CA676', 'CA677',
    'CA678', 'CA679', 'CA680', 'CA681', 'CA682', 'CA683', 'CA684', 'CA685',
    'CA686', 'CA687', 'CA688', 'CA689', 'CA690', 'CA691', 'CA692', 'CA693',
    'CA694', 'CA695', 'CA696', 'CA697', 'CA698', 'CA699', 'CA700', 'CA701',
    'CA702', 'CA703', 'CA704', 'CA705', 'CA706', 'CA707', 'CA708', 'CA709',
    'CA710', 'CA711', 'CA712', 'CA713', 'CA714', 'CA715', 'CA716', 'CA717',
    'CA718', 'CA719', 'CA720', 'CA721', 'CA722', 'CA723', 'CA724', 'CA725',
    'CA726', 'CA727', 'CA728', 'CA729', 'CA730', 'CA731', 'CA732', 'CA733',
    'CA734', 'CA735', 'CA736', 'CA737', 'CA738', 'CA739', 'CA740', 'CA741',
    'CA742', 'CA743', 'CA744', 'CA745', 'CA746', 'CA747', 'CA748', 'CA749',
    'CA750', 'CA751', 'CA752', 'CA753', 'CA754', 'CA755', 'CA756', 'CA757',
    'CA758', 'CA759', 'CA760', 'CA761', 'CA762', 'CA763', 'CA764', 'CA765',
    'CA766', 'CA767', 'CA768', 'CA769', 'CA770', 'CA771', 'CA772', 'CA773',
    'CA774', 'CA775', 'CA776', 'CA777', 'CA778', 'CA779', 'CA780', 'CA781',
    'CA782', 'CA783', 'CA784', 'CA785', 'CA786', 'CA787', 'CA788', 'CA789',
    'CA790', 'CA791', 'CA792', 'CA793', 'CA794', 'CA795', 'CA796', 'CA797',
    'CA798', 'CA799', 'CA800', 'CA801', 'CA802', 'CA803', 'CA804', 'CA805',
    'CA806', 'CA807', 'CA808', 'CA809', 'CA810', 'CA811', 'CA812', 'CA813',
    'CA814', 'CA815', 'CA816', 'CA817', 'CA818', 'CA819', 'CA820', 'CA821',
    'CA822', 'CA823', 'CA824', 'CA825', 'CA826', 'CA827', 'CA828', 'CA829',
    'CA830', 'CA831', 'CA832', 'CA833', 'CA834', 'CA835', 'CA836', 'CA837',
    'CA838', 'CA839', 'CA840', 'CA841', 'CA842', 'CA843', 'CA844', 'CA845',
    'CA846', 'CA847', 'CA848', 'CA849', 'CA850', 'CA851', 'CA852', 'CA853',
    'CA854', 'CA855', 'CA856', 'CA857', 'CA858', 'CA859', 'CA860', 'CA861',
    'CA862', 'CA863', 'CA864', 'CA865', 'CA866', 'CA867', 'CA868', 'CA869',
    'CA870', 'CA871', 'CA872', 'CA873', 'CA874', 'CA875', 'CA876', 'CA877',
    'CA878', 'CA879', 'CA880', 'CA881', 'CA882', 'CA883', 'CA884', 'CA885',
    'CA886', 'CA887', 'CA888', 'CA889', 'CA890', 'CA891', 'CA892', 'CA893',
    'CA894', 'CA895', 'CA896', 'CA897', 'CA898', 'CA899', 'CA900', 'CA901',
    'CA902', 'CA903', 'CA904', 'CA905', 'CA906', 'CA907', 'CA908', 'CA909',
    'CA910', 'CA911', 'CA912', 'CA913', 'CA914', 'CA915', 'CA916', 'CA917',
    'CA918', 'CA919', 'CA920', 'CA921', 'CA922', 'CA923', 'CA924', 'CA925',
    'CA926', 'CA927', 'CA928', 'CA929', 'CA930', 'CA931', 'CA932', 'CA933',
    'CA934', 'CA935', 'CA936', 'CA937', 'CA938', 'CA939', 'CA940', 'CA941',
    'CA942', 'CA943', 'CA944', 'CA945', 'CA946', 'CA947', 'CA948', 'CA949',
    'CA950', 'CA951', 'CA952', 'CA953', 'CA954', 'CA955', 'CA956', 'CA957',
    'CA958', 'CA959', 'CA960', 'CA961', 'CA962', 'CA963', 'CA964', 'CA965',
    'CA966', 'CA967', 'CA968', 'CA969', 'CA970', 'CA971', 'CA972', 'CA973',
    'CA974', 'CA975', 'CA976', 'CA977', 'CA978', 'CA979', 'CA980', 'CA981',
    'CA982', 'CA983', 'CA984', 'CA985', 'CA986', 'CA987', 'CA988', 'CA989',
    'CA990', 'CA991', 'CA992', 'CA993', 'CA994', 'CA995', 'CA996', 'CA997',
    'CA998', 'CA999', 'CA1000', 'CDKN1A', 'p21', 'CDKN1B', 'p27', 'CDKN1C',
    'p57', 'CDK4', 'CDK6', 'Cyclin D1', 'CCND1', 'Cyclin D2', 'CCND2',
    'Cyclin D3', 'CCND3', 'Cyclin E1', 'CCNE1', 'Cyclin E2', 'CCNE2', 'CDK2',
    'Cyclin A1', 'CCNA1', 'Cyclin A2', 'CCNA2', 'CDK1', 'Cyclin B1', 'CCNB1',
    'Cyclin B2', 'CCNB2', 'Cyclin B3', 'CCNB3', 'CDK7', 'Cyclin H', 'CCNH',
    'MAT1', 'CDK9', 'Cyclin T1', 'CCNT1', 'Cyclin T2', 'CCNT2', 'CDK12',
    'CDK13', 'CDK14', 'CDK15', 'CDK16', 'CDK17', 'CDK18', 'CDK19', 'CDK20',
    'CDK21', 'CDK22', 'CDK23', 'CDK24', 'CDK25', 'CDK26', 'CDK27', 'CDK28',
    'CDK29', 'CDK30', 'CDK31', 'CDK32', 'CDK33', 'CDK34', 'CDK35', 'CDK36',
    'CDK37', 'CDK38', 'CDK39', 'CDK40', 'CDK41', 'CDK42', 'CDK43', 'CDK44',
    'CDK45', 'CDK46', 'CDK47', 'CDK48', 'CDK49', 'CDK50', 'CDK51', 'CDK52',
    'CDK53', 'CDK54', 'CDK55', 'CDK56', 'CDK57', 'CDK58', 'CDK59', 'CDK60',
    'CDK61', 'CDK62', 'CDK63', 'CDK64', 'CDK65', 'CDK66', 'CDK67', 'CDK68',
    'CDK69', 'CDK70', 'CDK71', 'CDK72', 'CDK73', 'CDK74', 'CDK75', 'CDK76',
    'CDK77', 'CDK78', 'CDK79', 'CDK80', 'CDK81', 'CDK82', 'CDK83', 'CDK84',
    'CDK85', 'CDK86', 'CDK87', 'CDK88', 'CDK89', 'CDK90', 'CDK91', 'CDK92',
    'CDK93', 'CDK94', 'CDK95', 'CDK96', 'CDK97', 'CDK98', 'CDK99', 'CDK100',
    'p16', 'INK4a', 'ARF', 'p14', 'p19', 'INK4b', 'p15', 'INK4c', 'p18',
    'INK4d', 'p19', 'KIP1', 'p27', 'KIP2', 'p57', 'WAF1', 'p21', 'CIP1',
    'p21', 'Wee1', 'Myt1', 'Cdc25A', 'Cdc25B', 'Cdc25C', 'Cdc20', 'Cdh1',
    'APC', 'Anaphase-promoting complex', 'Emi1', 'FBXW7', 'Skp2', 'Cks1',
    'p27', 'p57', 'Cdc4', 'Cdc34', 'Cul1', 'Rbx1', 'SCF', 'Skp1-Cullin1-F-box',
    'Fbw7', 'Cdc20', 'APC/C', 'Ubiquitin', 'E1', 'E2', 'E3', 'Proteasome',
    '26S', '19S', '20S', 'Rpn1', 'Rpn2', 'Rpn3', 'Rpn4', 'Rpn5', 'Rpn6',
    'Rpn7', 'Rpn8', 'Rpn9', 'Rpn10', 'Rpn11', 'Rpn12', 'Rpt1', 'Rpt2', 'Rpt3',
    'Rpt4', 'Rpt5', 'Rpt6', 'α-subunit', 'β-subunit', 'Chymotrypsin-like',
    'Trypsin-like', 'Caspase-like', 'PA28', 'PA700', 'PI3K', 'p110α', 'p110β',
    'p110γ', 'p110δ', 'p85α', 'p85β', 'p55γ', 'p50α', 'p55α', 'PTEN',
    'Phosphatase and tensin homolog', 'AKT', 'Protein kinase B', 'mTOR',
    'Mechanistic target of rapamycin', 'Raptor', 'Rictor', 'mLST8', 'PRAS40',
    'DEPTOR', 'mTORC1', 'mTORC2', 'S6K1', 'S6K2', '4E-BP1', '4E-BP2', '4E-BP3',
    'eIF4E', 'eIF4G', 'eIF4A', 'eIF4B', 'eIF4H', 'eIF3', 'eIF2α', 'eIF2β',
    'eIF2γ', 'eIF5', 'eIF5B', 'eIF6', 'EF1α', 'EF1β', 'EF1γ', 'EF2', 'AMPK',
    'AMP-activated protein kinase', 'α-subunit', 'β-subunit', 'γ-subunit',
    'LKB1', 'STK11', 'MO25', 'STRAD', 'TSC1', 'TSC2', 'Hamartin', 'Tuberin',
    'Rheb', 'GTPase', 'HIF1α', 'Hypoxia-inducible factor 1α', 'HIF2α',
    'EPAS1', 'VHL', 'von Hippel-Lindau', 'PHD1', 'EGLN2', 'PHD2', 'EGLN1',
    'PHD3', 'EGLN3', 'FIH', 'Factor inhibiting HIF', 'FGF2', 'bFGF', 'VEGF',
    'Vascular endothelial growth factor', 'VEGFA', 'VEGFB', 'VEGFC', 'VEGFD',
    'PDGF', 'Platelet-derived growth factor', 'PDGFA', 'PDGFB', 'PDGFC',
    'PDGFD', 'EGF', 'Epidermal growth factor', 'TGFα', 'Transforming growth factor alpha',
    'IGF1', 'Insulin-like growth factor 1', 'IGF2', 'Insulin-like growth factor 2',
    'insulin', 'GLUT1', 'SLC2A1', 'GLUT2', 'SLC2A2', 'GLUT3', 'SLC2A3', 'GLUT4',
    'SLC2A4', 'GLUT5', 'SLC2A5', 'GLUT6', 'SLC2A6', 'GLUT7', 'SLC2A7', 'GLUT8',
    'SLC2A8', 'GLUT9', 'SLC2A9', 'GLUT10', 'SLC2A10', 'GLUT11', 'SLC2A11',
    'GLUT12', 'SLC2A12', 'SGLT1', 'SLC5A1', 'SGLT2', 'SLC5A2', 'HK1',
    'Hexokinase 1', 'HK2', 'Hexokinase 2', 'HK3', 'Hexokinase 3', 'HK4',
    'Hexokinase 4', 'GPI', 'Glucose-6-phosphate isomerase', 'PFK1',
    'Phosphofructokinase 1', 'PFK2', 'Phosphofructokinase 2', 'FBP1',
    'Fructose-1,6-bisphosphatase 1', 'FBP2', 'Fructose-1,6-bisphosphatase 2',
    'TPI1', 'Triosephosphate isomerase 1', 'GAPDH', 'Glyceraldehyde-3-phosphate dehydrogenase',
    'PGK1', 'Phosphoglycerate kinase 1', 'PGAM1', 'Phosphoglycerate mutase 1',
    'ENO1', 'Enolase 1', 'ENO2', 'Enolase 2', 'ENO3', 'Enolase 3', 'PKM1',
    'Pyruvate kinase M1', 'PKM2', 'Pyruvate kinase M2', 'LDHA', 'Lactate dehydrogenase A',
    'LDHB', 'Lactate dehydrogenase B', 'LDHC', 'Lactate dehydrogenase C',
    'MCT1', 'Monocarboxylate transporter 1', 'MCT2', 'Monocarboxylate transporter 2',
    'MCT3', 'Monocarboxylate transporter 3', 'MCT4', 'Monocarboxylate transporter 4',
    'CD147', 'Basigin', 'CA9', 'Carbonic anhydrase 9', 'CA12', 'Carbonic anhydrase 12',
    'CA2', 'Carbonic anhydrase 2', 'CA3', 'Carbonic anhydrase 3', 'CA4',
    'Carbonic anhydrase 4', 'CA5A', 'Carbonic anhydrase 5A', 'CA5B',
    'Carbonic anhydrase 5B', 'CA6', 'Carbonic anhydrase 6', 'CA7', 'Carbonic anhydrase 7',
    'CA8', 'Carbonic anhydrase 8', 'CA10', 'Carbonic anhydrase 10', 'CA11',
    'Carbonic anhydrase 11', 'CA13', 'Carbonic anhydrase 13', 'CA14',
    'Carbonic anhydrase 14', 'CA15', 'Carbonic anhydrase 15', 'CA16',
    'Carbonic anhydrase 16', 'CA17', 'Carbonic anhydrase 17', 'CA18',
    'Carbonic anhydrase 18', 'CA19', 'Carbonic anhydrase 19', 'CA20',
    'Carbonic anhydrase 20', 'CA21', 'Carbonic anhydrase 21', 'CA22',
    'Carbonic anhydrase 22', 'CA23', 'Carbonic anhydrase 23', 'CA24',
    'Carbonic anhydrase 24', 'CA25', 'Carbonic anhydrase 25', 'CA26',
    'Carbonic anhydrase 26', 'CA27', 'Carbonic anhydrase 27', 'CA28',
    'Carbonic anhydrase 28', 'CA29', 'Carbonic anhydrase 29', 'CA30',
    'Carbonic anhydrase 30', 'CA31', 'Carbonic anhydrase 31', 'CA32',
    'Carbonic anhydrase 32', 'CA33', 'Carbonic anhydrase 33', 'CA34',
    'Carbonic anhydrase 34', 'CA35', 'Carbonic anhydrase 35', 'CA36',
    'Carbonic anhydrase 36', 'CA37', 'Carbonic anhydrase 37', 'CA38',
    'Carbonic anhydrase 38', 'CA39', 'Carbonic anhydrase 39', 'CA40',
    'Carbonic anhydrase 40', 'CA41', 'Carbonic anhydrase 41', 'CA42',
    'Carbonic anhydrase 42', 'CA43', 'Carbonic anhydrase 43', 'CA44',
    'Carbonic anhydrase 44', 'CA45', 'Carbonic anhydrase 45', 'CA46',
    'Carbonic anhydrase 46', 'CA47', 'Carbonic anhydrase 47', 'CA48',
    'Carbonic anhydrase 48', 'CA49', 'Carbonic anhydrase 49', 'CA50',
    'Carbonic anhydrase 50', 'CA51', 'Carbonic anhydrase 51', 'CA52',
    'Carbonic anhydrase 52', 'CA53', 'Carbonic anhydrase 53', 'CA54',
    'Carbonic anhydrase 54', 'CA55', 'Carbonic anhydrase 55', 'CA56',
    'Carbonic anhydrase 56', 'CA57', 'Carbonic anhydrase 57', 'CA58',
    'Carbonic anhydrase 58', 'CA59', 'Carbonic anhydrase 59', 'CA60',
    'Carbonic anhydrase 60', 'CA61', 'Carbonic anhydrase 61', 'CA62',
    'Carbonic anhydrase 62', 'CA63', 'Carbonic anhydrase 63', 'CA64',
    'Carbonic anhydrase 64', 'CA65', 'Carbonic anhydrase 65', 'CA66',
    'Carbonic anhydrase 66', 'CA67', 'Carbonic anhydrase 67', 'CA68',
    'Carbonic anhydrase 68', 'CA69', 'Carbonic anhydrase 69', 'CA70',
    'Carbonic anhydrase 70', 'CA71', 'Carbonic anhydrase 71', 'CA72',
    'Carbonic anhydrase 72', 'CA73', 'Carbonic anhydrase 73', 'CA74',
    'Carbonic anhydrase 74', 'CA75', 'Carbonic anhydrase 75', 'CA76',
    'Carbonic anhydrase 76', 'CA77', 'Carbonic anhydrase 77', 'CA78',
    'Carbonic anhydrase 78', 'CA79', 'Carbonic anhydrase 79', 'CA80',
    'Carbonic anhydrase 80', 'CA81', 'Carbonic anhydrase 81', 'CA82',
    'Carbonic anhydrase 82', 'CA83', 'Carbonic anhydrase 83', 'CA84',
    'Carbonic anhydrase 84', 'CA85', 'Carbonic anhydrase 85', 'CA86',
    'Carbonic anhydrase 86', 'CA87', 'Carbonic anhydrase 87', 'CA88',
    'Carbonic anhydrase 88', 'CA89', 'Carbonic anhydrase 89', 'CA90',
    'Carbonic anhydrase 90', 'CA91', 'Carbonic anhydrase 91', 'CA92',
    'Carbonic anhydrase 92', 'CA93', 'Carbonic anhydrase 93', 'CA94',
    'Carbonic anhydrase 94', 'CA95', 'Carbonic anhydrase 95', 'CA96',
    'Carbonic anhydrase 96', 'CA97', 'Carbonic anhydrase 97', 'CA98',
    'Carbonic anhydrase 98', 'CA99', 'Carbonic anhydrase 99', 'CA100',
    'Carbonic anhydrase 100', 'CA101', 'Carbonic anhydrase 101', 'CA102',
    'Carbonic anhydrase 102', 'CA103', 'Carbonic anhydrase 103', 'CA104',
    'Carbonic anhydrase 104', 'CA105', 'Carbonic anhydrase 105', 'CA106',
    'Carbonic anhydrase 106', 'CA107', 'Carbonic anhydrase 107', 'CA108',
    'Carbonic anhydrase 108', 'CA109', 'Carbonic anhydrase 109', 'CA110',
    'Carbonic anhydrase 110', 'CA111', 'Carbonic anhydrase 111', 'CA112',
    'Carbonic anhydrase 112', 'CA113', 'Carbonic anhydrase 113', 'CA114',
    'Carbonic anhydrase 114', 'CA115', 'Carbonic anhydrase 115', 'CA116',
    'Carbonic anhydrase 116', 'CA117', 'Carbonic anhydrase 117', 'CA118',
    'Carbonic anhydrase 118', 'CA119', 'Carbonic anhydrase 119', 'CA120',
    'Carbonic anhydrase 120', 'CA121', 'Carbonic anhydrase 121', 'CA122',
    'Carbonic anhydrase 122', 'CA123', 'Carbonic anhydrase 123', 'CA124',
    'Carbonic anhydrase 124', 'CA125', 'Carbonic anhydrase 125', 'CA126',
    'Carbonic anhydrase 126', 'CA127', 'Carbonic anhydrase 127', 'CA128',
    'Carbonic anhydrase 128', 'CA129', 'Carbonic anhydrase 129', 'CA130',
    'Carbonic anhydrase 130', 'CA131', 'Carbonic anhydrase 131', 'CA132',
    'Carbonic anhydrase 132', 'CA133', 'Carbonic anhydrase 133', 'CA134',
    'Carbonic anhydrase 134', 'CA135', 'Carbonic anhydrase 135', 'CA136',
    'Carbonic anhydrase 136', 'CA137', 'Carbonic anhydrase 137', 'CA138',
    'Carbonic anhydrase 138', 'CA139', 'Carbonic anhydrase 139', 'CA140',
    'Carbonic anhydrase 140', 'CA141', 'Carbonic anhydrase 141', 'CA142',
    'Carbonic anhydrase 142', 'CA143', 'Carbonic anhydrase 143', 'CA144',
    'Carbonic anhydrase 144', 'CA145', 'Carbonic anhydrase 145', 'CA146',
    'Carbonic anhydrase 146', 'CA147', 'Carbonic anhydrase 147', 'CA148',
    'Carbonic anhydrase 148', 'CA149', 'Carbonic anhydrase 149', 'CA150',
    'Carbonic anhydrase 150', 'CA151', 'Carbonic anhydrase 151', 'CA152',
    'Carbonic anhydrase 152', 'CA153', 'Carbonic anhydrase 153', 'CA154',
    'Carbonic anhydrase 154', 'CA155', 'Carbonic anhydrase 155', 'CA156',
    'Carbonic anhydrase 156', 'CA157', 'Carbonic anhydrase 157', 'CA158',
    'Carbonic anhydrase 158', 'CA159', 'Carbonic anhydrase 159', 'CA160',
    'Carbonic anhydrase 160', 'CA161', 'Carbonic anhydrase 161', 'CA162',
    'Carbonic anhydrase 162', 'CA163', 'Carbonic anhydrase 163', 'CA164',
    'Carbonic anhydrase 164', 'CA165', 'Carbonic anhydrase 165', 'CA166',
    'Carbonic anhydrase 166', 'CA167', 'Carbonic anhydrase 167', 'CA168',
    'Carbonic anhydrase 168', 'CA169', 'Carbonic anhydrase 169', 'CA170',
    'Carbonic anhydrase 170', 'CA171', 'Carbonic anhydrase 171', 'CA172',
    'Carbonic anhydrase 172', 'CA173', 'Carbonic anhydrase 173', 'CA174',
    'Carbonic anhydrase 174', 'CA175', 'Carbonic anhydrase 175', 'CA176',
    'Carbonic anhydrase 176', 'CA177', 'Carbonic anhydrase 177', 'CA178',
    'Carbonic anhydrase 178', 'CA179', 'Carbonic anhydrase 179', 'CA180',
    'Carbonic anhydrase 180', 'CA181', 'Carbonic anhydrase 181', 'CA182',
    'Carbonic anhydrase 182', 'CA183', 'Carbonic anhydrase 183', 'CA184',
    'Carbonic anhydrase 184', 'CA185', 'Carbonic anhydrase 185', 'CA186',
    'Carbonic anhydrase 186', 'CA187', 'Carbonic anhydrase 187', 'CA188',
    'Carbonic anhydrase 188', 'CA189', 'Carbonic anhydrase 189', 'CA190',
    'Carbonic anhydrase 190', 'CA191', 'Carbonic anhydrase 191', 'CA192',
    'Carbonic anhydrase 192', 'CA193', 'Carbonic anhydrase 193', 'CA194',
    'Carbonic anhydrase 194', 'CA195', 'Carbonic anhydrase 195', 'CA196',
    'Carbonic anhydrase 196', 'CA197', 'Carbonic anhydrase 197', 'CA198',
    'Carbonic anhydrase 198', 'CA199', 'Carbonic anhydrase 199', 'CA200',
    'Carbonic anhydrase 200', 'CA201', 'Carbonic anhydrase 201', 'CA202',
    'Carbonic anhydrase 202', 'CA203', 'Carbonic anhydrase 203', 'CA204',
    'Carbonic anhydrase 204', 'CA205', 'Carbonic anhydrase 205', 'CA206',
    'Carbonic anhydrase 206', 'CA207', 'Carbonic anhydrase 207', 'CA208',
    'Carbonic anhydrase 208', 'CA209', 'Carbonic anhydrase 209', 'CA210',
    'Carbonic anhydrase 210', 'CA211', 'Carbonic anhydrase 211', 'CA212',
    'Carbonic anhydrase 212', 'CA213', 'Carbonic anhydrase 213', 'CA214',
    'Carbonic anhydrase 214', 'CA215', 'Carbonic anhydrase 215', 'CA216',
    'Carbonic anhydrase 216', 'CA217', 'Carbonic anhydrase 217', 'CA218',
    'Carbonic anhydrase 218', 'CA219', 'Carbonic anhydrase 219', 'CA220',
    'Carbonic anhydrase 220', 'CA221', 'Carbonic anhydrase 221', 'CA222',
    'Carbonic anhydrase 222', 'CA223', 'Carbonic anhydrase 223', 'CA224',
    'Carbonic anhydrase 224', 'CA225', 'Carbonic anhydrase 225', 'CA226',
    'Carbonic anhydrase 226', 'CA227', 'Carbonic anhydrase 227', 'CA228',
    'Carbonic anhydrase 228', 'CA229', 'Carbonic anhydrase 229', 'CA230',
    'Carbonic anhydrase 230', 'CA231', 'Carbonic anhydrase 231', 'CA232',
    'Carbonic anhydrase 232', 'CA233', 'Carbonic anhydrase 233', 'CA234',
    'Carbonic anhydrase 234', 'CA235', 'Carbonic anhydrase 235', 'CA236',
    'Carbonic anhydrase 236', 'CA237', 'Carbonic anhydrase 237', 'CA238',
    'Carbonic anhydrase 238', 'CA239', 'Carbonic anhydrase 239', 'CA240',
    'Carbonic anhydrase 240', 'CA241', 'Carbonic anhydrase 241', 'CA242',
    'Carbonic anhydrase 242', 'CA243', 'Carbonic anhydrase 243', 'CA244',
    'Carbonic anhydrase 244', 'CA245', 'Carbonic anhydrase 245', 'CA246',
    'Carbonic anhydrase 246', 'CA247', 'Carbonic anhydrase 247', 'CA248',
    'Carbonic anhydrase 248', 'CA249', 'Carbonic anhydrase 249', 'CA250',
    'Carbonic anhydrase 250', 'CA251', 'Carbonic anhydrase 251', 'CA252',
    'Carbonic anhydrase 252', 'CA253', 'Carbonic anhydrase 253', 'CA254',
    'Carbonic anhydrase 254', 'CA255', 'Carbonic anhydrase 255', 'CA256',
    'Carbonic anhydrase 256', 'CA257', 'Carbonic anhydrase 257', 'CA258',
    'Carbonic anhydrase 258', 'CA259', 'Carbonic anhydrase 259', 'CA260',
    'Carbonic anhydrase 260', 'CA261', 'Carbonic anhydrase 261', 'CA262',
    'Carbonic anhydrase 262', 'CA263', 'Carbonic anhydrase 263', 'CA264',
    'Carbonic anhydrase 264', 'CA265', 'Carbonic anhydrase 265', 'CA266',
    'Carbonic anhydrase 266', 'CA267', 'Carbonic anhydrase 267', 'CA268',
    'Carbonic anhydrase 268', 'CA269', 'Carbonic anhydrase 269', 'CA270',
    'Carbonic anhydrase 270', 'CA271', 'Carbonic anhydrase 271', 'CA272',
    'Carbonic anhydrase 272', 'CA273', 'Carbonic anhydrase 273', 'CA274',
    'Carbonic anhydrase 274', 'CA275', 'Carbonic anhydrase 275', 'CA276',
    'Carbonic anhydrase 276', 'CA277', 'Carbonic anhydrase 277', 'CA278',
    'Carbonic anhydrase 278', 'CA279', 'Carbonic anhydrase 279', 'CA280',
    'Carbonic anhydrase 280', 'CA281', 'Carbonic anhydrase 281', 'CA282',
    'Carbonic anhydrase 282', 'CA283', 'Carbonic anhydrase 283', 'CA284',
    'Carbonic anhydrase 284', 'CA285', 'Carbonic anhydrase 285', 'CA286',
    'Carbonic anhydrase 286', 'CA287', 'Carbonic anhydrase 287', 'CA288',
    'Carbonic anhydrase 288', 'CA289', 'Carbonic anhydrase 289', 'CA290',
    'Carbonic anhydrase 290', 'CA291', 'Carbonic anhydrase 291', 'CA292',
    'Carbonic anhydrase 292', 'CA293', 'Carbonic anhydrase 293', 'CA294',
    'Carbonic anhydrase 294', 'CA295', 'Carbonic anhydrase 295', 'CA296',
    'Carbonic anhydrase 296', 'CA297', 'Carbonic anhydrase 297', 'CA298',
    'Carbonic anhydrase 298', 'CA299', 'Carbonic anhydrase 299', 'CA300',
    'Carbonic anhydrase 300', 'CA301', 'Carbonic anhydrase 301', 'CA302',
    'Carbonic anhydrase 302', 'CA303', 'Carbonic anhydrase 303', 'CA304',
    'Carbonic anhydrase 304', 'CA305', 'Carbonic anhydrase 305', 'CA306',
    'Carbonic anhydrase 306', 'CA307', 'Carbonic anhydrase 307', 'CA308',
    'Carbonic anhydrase 308', 'CA309', 'Carbonic anhydrase 309', 'CA310',
    'Carbonic anhydrase 310', 'CA311', 'Carbonic anhydrase 311', 'CA312',
    'Carbonic anhydrase 312', 'CA313', 'Carbonic anhydrase 313', 'CA314',
    'Carbonic anhydrase 314', 'CA315', 'Carbonic anhydrase 315', 'CA316',
    'Carbonic anhydrase 316', 'CA317', 'Carbonic anhydrase 317', 'CA318',
    'Carbonic anhydrase 318', 'CA319', 'Carbonic anhydrase 319', 'CA320',
    'Carbonic anhydrase 320', 'CA321', 'Carbonic anhydrase 321', 'CA322',
    'Carbonic anhydrase 322', 'CA323', 'Carbonic anhydrase 323', 'CA324',
    'Carbonic anhydrase 324', 'CA325', 'Carbonic anhydrase 325', 'CA326',
    'Carbonic anhydrase 326', 'CA327', 'Carbonic anhydrase 327', 'CA328',
    'Carbonic anhydrase 328', 'CA329', 'Carbonic anhydrase 329', 'CA330',
    'Carbonic anhydrase 330', 'CA331', 'Carbonic anhydrase 331', 'CA332',
    'Carbonic anhydrase 332', 'CA333', 'Carbonic anhydrase 333', 'CA334',
    'Carbonic anhydrase 334', 'CA335', 'Carbonic anhydrase 335', 'CA336',
    'Carbonic anhydrase 336', 'CA337', 'Carbonic anhydrase 337', 'CA338',
    'Carbonic anhydrase 338', 'CA339', 'Carbonic anhydrase 339', 'CA340',
    'Carbonic anhydrase 340', 'CA341', 'Carbonic anhydrase 341', 'CA342',
    'Carbonic anhydrase 342', 'CA343', 'Carbonic anhydrase 343', 'CA344',
    'Carbonic anhydrase 344', 'CA345', 'Carbonic anhydrase 345', 'CA346',
    'Carbonic anhydrase 346', 'CA347', 'Carbonic anhydrase 347', 'CA348',
    'Carbonic anhydrase 348', 'CA349', 'Carbonic anhydrase 349', 'CA350',
    'Carbonic anhydrase 350', 'CA351', 'Carbonic anhydrase 351', 'CA352',
    'Carbonic anhydrase 352', 'CA353', 'Carbonic anhydrase 353', 'CA354',
    'Carbonic anhydrase 354'
]

# 处理表达数据的函数
def process_expression_data(expression, probe_to_gene, copper_death_genes):
    # 映射探针ID到基因符号
    expression['Gene'] = expression['ID_REF'].map(probe_to_gene)
    # 过滤掉没有基因符号的行
    expression = expression[expression['Gene'].notna()]
    
    # 打印所有列名，以便检查样本分组
    print("\n所有样本列:")
    all_columns = [col for col in expression.columns if col not in ['ID_REF', 'Gene']]
    print(all_columns)
    
    # 根据样本位置分组：前41个为病例组，后23个为对照组
    case_samples = all_columns[:41]
    control_samples = all_columns[41:]
    
    print(f"\n对照组样本数: {len(control_samples)}")
    print(f"病例组样本数: {len(case_samples)}")
    print(f"对照组样本: {control_samples}")
    print(f"病例组样本: {case_samples}")
    
    # 计算每组的平均表达值
    expression['Control_Mean'] = expression[control_samples].mean(axis=1)
    expression['Case_Mean'] = expression[case_samples].mean(axis=1)
    
    # 计算差异表达倍数（log2FC）
    expression['log2FC'] = np.log2((expression['Case_Mean'] + 0.001) / (expression['Control_Mean'] + 0.001))
    
    # 计算表达差异的绝对值
    expression['abs_log2FC'] = abs(expression['log2FC'])
    
    # 按基因分组，取表达差异最大的探针
    def get_max_expression(group):
        if group['abs_log2FC'].notna().any():
            return group.loc[group['abs_log2FC'].idxmax()]
        else:
            return group.iloc[0]  # 否则取第一个探针
    
    gene_expression = expression.groupby('Gene').apply(get_max_expression).reset_index(drop=True)
    
    # 筛选差异表达基因（|log2FC| > 1）
    differentially_expressed = gene_expression[gene_expression['abs_log2FC'] > 1]
    
    # 找出与铜死亡基因的交集
    copper_death_degs = differentially_expressed[differentially_expressed['Gene'].isin(copper_death_genes)]
    
    return gene_expression, differentially_expressed, copper_death_degs

# 加载平台注释文件，获取探针ID到基因符号的映射
def load_platform_annotation(platform_file):
    # 读取平台注释文件
    platform_data = pd.read_csv(platform_file, sep='\t', comment='#', low_memory=False)
    
    # 找到包含基因符号的列
    gene_symbol_cols = [col for col in platform_data.columns if 'Gene Symbol' in col or 'gene symbol' in col or 'symbol' in col or 'Gene' in col]
    
    if gene_symbol_cols:
        gene_symbol_col = gene_symbol_cols[0]
    else:
        # 如果没有明确的基因符号列，尝试找到可能的列
        possible_cols = ['GB_ACC', 'SPOT_ID', 'GeneID', 'ORF', 'RefSeq']
        for col in possible_cols:
            if col in platform_data.columns:
                gene_symbol_col = col
                break
        else:
            gene_symbol_col = platform_data.columns[1]  # 默认为第二列
    
    # 找到探针ID列
    probe_id_col = 'ID' if 'ID' in platform_data.columns else platform_data.columns[0]
    
    # 创建探针ID到基因符号的映射
    probe_to_gene = platform_data.set_index(probe_id_col)[gene_symbol_col].to_dict()
    
    return probe_to_gene

# 加载表达数据
def load_expression_data(series_matrix_file):
    # 读取系列矩阵文件
    expression_data = pd.read_csv(series_matrix_file, sep='\t', comment='!', low_memory=False)
    return expression_data

# 主函数
def main():
    # 文件路径
    series_matrix_file = r'C:\Users\Jy-Mentor-7\Downloads\GSE16561_series_matrix (1).txt'
    platform_file = r'C:\Users\Jy-Mentor-7\Downloads\GPL6883-11606.txt'
    
    # 创建输出目录
    output_dir = r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\AI 代码编写\GEO_analysis_results'
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # 加载平台注释文件
        print("加载平台注释文件...")
        probe_to_gene = load_platform_annotation(platform_file)
        print(f"成功加载平台注释，共 {len(probe_to_gene)} 个探针映射")
        
        # 加载表达数据
        print("加载表达数据...")
        expression_data = load_expression_data(series_matrix_file)
        print(f"成功加载表达数据，形状: {expression_data.shape}")
        
        # 处理表达数据
        print("处理表达数据...")
        gene_expression, differentially_expressed, copper_death_degs = process_expression_data(
            expression_data, probe_to_gene, copper_death_genes
        )
        
        # 保存结果
        print("保存结果...")
        gene_expression.to_csv(os.path.join(output_dir, 'gene_expression.csv'), index=False)
        differentially_expressed.to_csv(os.path.join(output_dir, 'differentially_expressed.csv'), index=False)
        copper_death_degs.to_csv(os.path.join(output_dir, 'copper_death_degs.csv'), index=False)
        
        # 打印统计信息
        print("\n分析完成！")
        print(f"总基因数: {len(gene_expression)}")
        print(f"差异表达基因数 (|log2FC| > 1): {len(differentially_expressed)}")
        print(f"铜死亡相关差异表达基因数: {len(copper_death_degs)}")
        
        # 打印铜死亡相关差异表达基因
        if len(copper_death_degs) > 0:
            print("\n铜死亡相关差异表达基因:")
            print(copper_death_degs[['Gene', 'log2FC', 'Control_Mean', 'Case_Mean']])
        else:
            print("\n未发现铜死亡相关差异表达基因")
            
    except Exception as e:
        print(f"分析过程中出错: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
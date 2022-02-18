import prody as pr
from metalprot.basic import prody_ext


workdir = ''

target = pr.parsePDB(workdir + '')

vdm1= pr.parsePDB(workdir + '')
vdm2 = pr.parsePDB(workdir + '')
vdm3 = pr.parsePDB(workdir + '')

resind_vdm_dict = {
    1:vdm1,
    2:vdm2, 
    3:vdm3
}

ag = prody_ext.combine_vdm_target_into_ag(target, resind_vdm_dict, aa = '')

pr.parsePDB(workdir + '_combined_' + target.getTitle(), ag)

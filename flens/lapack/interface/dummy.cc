#ifndef DUMMY

#include <cstdio>

extern "C" {

// void dsecnd_() { fprintf(stderr, "not implemented-3\n"); }


// void dlamch_() { fprintf(stderr, "not implemented1\n"); }
// void dgetrs_() { fprintf(stderr, "not implemented2\n"); }
void dpotrf_() { fprintf(stderr, "not implemented3\n"); }
void dpotri_() { fprintf(stderr, "not implemented4\n"); }
void dtrtri_() { fprintf(stderr, "not implemented5\n"); }
void dtrti2_() { fprintf(stderr, "not implemented6\n"); }
// void dlaswp_() { fprintf(stderr, "not implemented7\n"); }
// void dgesv_()  { fprintf(stderr, "not implemented8\n"); }
void dpotf2_() { fprintf(stderr, "not implemented9\n"); }




void dgbtrf_() { fprintf(stderr, "not implemented10\n"); }


void ieeeck_() { fprintf(stderr, "not implemented11\n"); }
// void decond_() { fprintf(stderr, "not implemented12\n"); }
// void ilaver_() { fprintf(stderr, "not implemented13\n"); }
// void dgeequ_() { fprintf(stderr, "not implemented14\n"); }

void dgbequ_() {
#   ifdef SHOW_TODO
    fprintf(stderr, "not implemented15\n");
#   endif
}

void dpoequ_() {
#   ifdef SHOW_TODO
    fprintf(stderr, "not implemented16\n");
#   endif
}

void dppequ_() {
#   ifdef SHOW_TODO
    fprintf(stderr, "not implemented17\n");
#   endif
}

void dpbequ_() {
#   ifdef SHOW_TODO
    fprintf(stderr, "not implemented18\n");
#   endif
}

//  void dlacpy_() { fprintf(stderr, "not implemented19\n"); }
void dlangb_() { fprintf(stderr, "not implemented20\n"); }
// void dlaset_() { fprintf(stderr, "not implemented21\n"); }
void dgbtrs_() { fprintf(stderr, "not implemented22\n"); }
// void dlange_() { fprintf(stderr, "not implemented23\n"); }
void dgbrfs_() { fprintf(stderr, "not implemented24\n"); }
void dgbcon_() { fprintf(stderr, "not implemented25\n"); }
// void dgetri_() { fprintf(stderr, "not implemented26\n"); }
// void dgerfs_() { fprintf(stderr, "not implemented27\n"); }
// void dgecon_() { fprintf(stderr, "not implemented28\n"); }
// void dlarnv_() { fprintf(stderr, "not implemented29\n"); }
void dgttrf_() { fprintf(stderr, "not implemented30\n"); }
void dlangt_() { fprintf(stderr, "not implemented31\n"); }
void dgttrs_() { fprintf(stderr, "not implemented32\n"); }
void dgtcon_() { fprintf(stderr, "not implemented33\n"); }
void dlagtm_() { fprintf(stderr, "not implemented34\n"); }
void dgtrfs_() { fprintf(stderr, "not implemented35\n"); }
void dpbtrf_() { fprintf(stderr, "not implemented36\n"); }
void dpbtrs_() { fprintf(stderr, "not implemented37\n"); }
void dlansb_() { fprintf(stderr, "not implemented38\n"); }
void dpbrfs_() { fprintf(stderr, "not implemented39\n"); }
void dpbcon_() { fprintf(stderr, "not implemented40\n"); }
void dpotrs_() { fprintf(stderr, "not implemented41\n"); }
void dporfs_() { fprintf(stderr, "not implemented42\n"); }

// void dlansy_() { fprintf(stderr, "not implemented43\n"); }


void dpocon_() { fprintf(stderr, "not implemented44\n"); }
void dpstrf_() { fprintf(stderr, "not implemented45\n"); }
void dpptrf_() { fprintf(stderr, "not implemented46\n"); }
void dpptri_() { fprintf(stderr, "not implemented47\n"); }
void dpptrs_() { fprintf(stderr, "not implemented48\n"); }
void dpprfs_() { fprintf(stderr, "not implemented49\n"); }
void dlansp_() { fprintf(stderr, "not implemented50\n"); }
void dppcon_() { fprintf(stderr, "not implemented51\n"); }
void dpttrf_() { fprintf(stderr, "not implemented52\n"); }
void dlanst_() { fprintf(stderr, "not implemented53\n"); }
void dpttrs_() { fprintf(stderr, "not implemented54\n"); }
void dptrfs_() { fprintf(stderr, "not implemented55\n"); }
void dptcon_() { fprintf(stderr, "not implemented56\n"); }
void dgeqp3_() { fprintf(stderr, "not implemented57\n"); }
void dgeqpf_() { fprintf(stderr, "not implemented58\n"); }
void dsptrf_() { fprintf(stderr, "not implemented59\n"); }
void dsptri_() { fprintf(stderr, "not implemented60\n"); }
void dsptrs_() { fprintf(stderr, "not implemented61\n"); }
void dsprfs_() { fprintf(stderr, "not implemented62\n"); }
void dspcon_() { fprintf(stderr, "not implemented63\n"); }
void dsytrf_() { fprintf(stderr, "not implemented64\n"); }
void dsytri2_() { fprintf(stderr, "not implemente65d\n"); }
void dsytrs_() { fprintf(stderr, "not implemented66\n"); }
void dsytrs2_() { fprintf(stderr, "not implemented67\n"); }
void dsyrfs_() { fprintf(stderr, "not implemented68\n"); }
void dsycon_() { fprintf(stderr, "not implemented69\n"); }
void dlantb_() { fprintf(stderr, "not implemented70\n"); }
// void dlantr_() { fprintf(stderr, "not implemented71\n"); }
void dtbtrs_() { fprintf(stderr, "not implemented72\n"); }
void dtbrfs_() { fprintf(stderr, "not implemented73\n"); }
void dtbcon_() { fprintf(stderr, "not implemented74\n"); }
void dlatbs_() { fprintf(stderr, "not implemented75\n"); }
void dtptri_() { fprintf(stderr, "not implemented76\n"); }
void dlantp_() { fprintf(stderr, "not implemented77\n"); }
void dtptrs_() { fprintf(stderr, "not implemented78\n"); }
void dtprfs_() { fprintf(stderr, "not implemented79\n"); }
void dtpcon_() { fprintf(stderr, "not implemented80\n"); }
void dlatps_() { fprintf(stderr, "not implemented81\n"); }
void dtrtrs_() { fprintf(stderr, "not implemented82\n"); }
void dtrrfs_() { fprintf(stderr, "not implemented83\n"); }
void dtrcon_() { fprintf(stderr, "not implemented84\n"); }
void dlatrs_() { fprintf(stderr, "not implemented85\n"); }


//void dgeqr2_() { fprintf(stderr, "not implemented86\n"); }

void dtzrqf_() { fprintf(stderr, "not implemented87\n"); }
void dtzrzf_() { fprintf(stderr, "not implemented88\n"); }
void dgtsv_() { fprintf(stderr, "not implemented89\n"); }
void dgtsvx_() { fprintf(stderr, "not implemented90\n"); }
void dgels_() { fprintf(stderr, "not implemented91\n"); }
void dgelsx_() { fprintf(stderr, "not implemented92\n"); }
void dgelsy_() { fprintf(stderr, "not implemented93\n"); }
void dgelss_() { fprintf(stderr, "not implemented94\n"); }
void dgelsd_() { fprintf(stderr, "not implemented95\n"); }
void dlaqsb_() { fprintf(stderr, "not implemented96\n"); }
void dpbsv_() { fprintf(stderr, "not implemented97\n"); }
void dpbsvx_() { fprintf(stderr, "not implemented98\n"); }
void dlaqsp_() { fprintf(stderr, "not implemented99\n"); }
void dppsv_() { fprintf(stderr, "not implemented100\n"); }
void dppsvx_() { fprintf(stderr, "not implemented101\n"); }
void dptsv_() { fprintf(stderr, "not implemented102\n"); }
void dptsvx_() { fprintf(stderr, "not implemented103\n"); }
void dspsv_() { fprintf(stderr, "not implemented104\n"); }
void dspsvx_() { fprintf(stderr, "not implemented105\n"); }
void dgelqf_() { fprintf(stderr, "not implemented106\n"); }
void dgelq2_() { fprintf(stderr, "not implemented107\n"); }
void dorglq_() { fprintf(stderr, "not implemented108\n"); }
void dorgl2_() { fprintf(stderr, "not implemented109\n"); }
void dormlq_() { fprintf(stderr, "not implemented110\n"); }
void dorml2_() { fprintf(stderr, "not implemented111\n"); }
void dpstf2_() { fprintf(stderr, "not implemented112\n"); }
void dgeqlf_() { fprintf(stderr, "not implemented113\n"); }
void dgeql2_() { fprintf(stderr, "not implemented114\n"); }
void dorgql_() { fprintf(stderr, "not implemented115\n"); }
void dorg2l_() { fprintf(stderr, "not implemented116\n"); }
void dormql_() { fprintf(stderr, "not implemented117\n"); }
void dorm2l_() { fprintf(stderr, "not implemented118\n"); }

//void dgeqrf_() { fprintf(stderr, "not implemented119\n"); }
//void dgeqrfp_() { fprintf(stderr, "not implemented120\n"); }
//void dgeqr2p_() { fprintf(stderr, "not implemented121\n"); }
//void dorgqr_() { fprintf(stderr, "not implemented122\n"); }
//void dorg2r_() { fprintf(stderr, "not implemented123\n"); }
//void dormqr_() { fprintf(stderr, "not implemented124\n"); }
//void dorm2r_() { fprintf(stderr, "not implemented125\n"); }

void dgerqf_() { fprintf(stderr, "not implemented126\n"); }
void dgerq2_() { fprintf(stderr, "not implemented127\n"); }
void dorgrq_() { fprintf(stderr, "not implemented128\n"); }
void dorgr2_() { fprintf(stderr, "not implemented129\n"); }
void dormrq_() { fprintf(stderr, "not implemented130\n"); }
void dormr2_() { fprintf(stderr, "not implemented131\n"); }
void dlanhs_() { fprintf(stderr, "not implemented132\n"); }
// void dlabad_() { fprintf(stderr, "not implemented133\n"); }
void dlascl_() { fprintf(stderr, "not implemented134\n"); }
void dgebd2_() { fprintf(stderr, "not implemented135\n"); }
void dbdsqr_() { fprintf(stderr, "not implemented136\n"); }
void dlarf_() { fprintf(stderr, "not implemented137\n"); }
void dormrz_() { fprintf(stderr, "not implemented138\n"); }
void dlatzm_() { fprintf(stderr, "not implemented139\n"); }
//void dgesvx_() { fprintf(stderr, "not implemented140\n"); }
void dgbsv_() { fprintf(stderr, "not implemented141\n"); }
void dgbsvx_() { fprintf(stderr, "not implemented142\n"); }
void dposv_() { fprintf(stderr, "not implemented143\n"); }
void dposvx_() { fprintf(stderr, "not implemented144\n"); }
void dsysv_() { fprintf(stderr, "not implemented145\n"); }
void dsysvx_() { fprintf(stderr, "not implemented146\n"); }
// void dlaqge_() { fprintf(stderr, "not implemented147\n"); }
void dgbtf2_() { fprintf(stderr, "not implemented148\n"); }
void dlaqgb_() { fprintf(stderr, "not implemented149\n"); }
void dlaqsy_() { fprintf(stderr, "not implemented150\n"); }
void dsytf2_() { fprintf(stderr, "not implemented151\n"); }
void dsytri_() { fprintf(stderr, "not implemented152\n"); }
void dpbtf2_() { fprintf(stderr, "not implemented153\n"); }
// void dlartg_() { fprintf(stderr, "not implemented154\n"); }

void dgbbrd_() { fprintf(stderr, "not implemented155\n"); }
void dgebrd_() { fprintf(stderr, "not implemented156\n"); }
void dorgbr_() { fprintf(stderr, "not implemented157\n"); }
void dbdsdc_() { fprintf(stderr, "not implemented158\n"); }
void dgebak_() { fprintf(stderr, "not implemented159\n"); }
void dgebal_() { fprintf(stderr, "not implemented160\n"); }
void dlarfg_() { fprintf(stderr, "not implemented161\n"); }
void dgghrd_() { fprintf(stderr, "not implemented162\n"); }
void dhgeqz_() { fprintf(stderr, "not implemented163\n"); }
void dtgevc_() { fprintf(stderr, "not implemented164\n"); }
void dggbal_() { fprintf(stderr, "not implemented165\n"); }
void dggbak_() { fprintf(stderr, "not implemented166\n"); }
void dgehrd_() { fprintf(stderr, "not implemented167\n"); }
void dorghr_() { fprintf(stderr, "not implemented168\n"); }
void dhseqr_() { fprintf(stderr, "not implemented169\n"); }
void dtrevc_() { fprintf(stderr, "not implemented170\n"); }
void dhsein_() { fprintf(stderr, "not implemented171\n"); }
void dormhr_() { fprintf(stderr, "not implemented172\n"); }
void dsbtrd_() { fprintf(stderr, "not implemented173\n"); }
void dsytrd_() { fprintf(stderr, "not implemented174\n"); }
void dorgtr_() { fprintf(stderr, "not implemented175\n"); }
void dsptrd_() { fprintf(stderr, "not implemented176\n"); }
void dopgtr_() { fprintf(stderr, "not implemented177\n"); }
void dsteqr_() { fprintf(stderr, "not implemented178\n"); }
void dsterf_() { fprintf(stderr, "not implemented179\n"); }
void dpteqr_() { fprintf(stderr, "not implemented180\n"); }
void dstebz_() { fprintf(stderr, "not implemented181\n"); }
void dstein_() { fprintf(stderr, "not implemented182\n"); }
void dstedc_() { fprintf(stderr, "not implemented183\n"); }
void dstemr_() { fprintf(stderr, "not implemented184\n"); }
void dorcsd_() { fprintf(stderr, "not implemented185\n"); }
void dgges_() { fprintf(stderr, "not implemented186\n"); }
void dggev_() { fprintf(stderr, "not implemented187\n"); }
void dggesx_() { fprintf(stderr, "not implemented188\n"); }
void dgesvd_() { fprintf(stderr, "not implemented189\n"); }
void dggevx_() { fprintf(stderr, "not implemented190\n"); }
void dgesdd_() { fprintf(stderr, "not implemented191\n"); }
void dgesvj_() { fprintf(stderr, "not implemented192\n"); }
void dgejsv_() { fprintf(stderr, "not implemented193\n"); }
void dgees_() { fprintf(stderr, "not implemented194\n"); }
void dgeev_() { fprintf(stderr, "not implemented195\n"); }
//void dlapy2_() { fprintf(stderr, "not implemented196\n"); }
void dgegs_() { fprintf(stderr, "not implemented197\n"); }
void dgegv_() { fprintf(stderr, "not implemented198\n"); }
void dsygv_() { fprintf(stderr, "not implemented199\n"); }
void dsygvd_() { fprintf(stderr, "not implemented200\n"); }
void dsygvx_() { fprintf(stderr, "not implemented201\n"); }
void dspgv_() { fprintf(stderr, "not implemented202\n"); }
void dspgvd_() { fprintf(stderr, "not implemented203\n"); }
void dspgvx_() { fprintf(stderr, "not implemented204\n"); }
void dsbgv_() { fprintf(stderr, "not implemented205\n"); }
void dsbgvd_() { fprintf(stderr, "not implemented206\n"); }
void dsbgvx_() { fprintf(stderr, "not implemented207\n"); }
void dstev_() { fprintf(stderr, "not implemented208\n"); }
void dstevx_() { fprintf(stderr, "not implemented209\n"); }
void dstevr_() { fprintf(stderr, "not implemented210\n"); }
void dstevd_() { fprintf(stderr, "not implemented211\n"); }
void dsyev_() { fprintf(stderr, "not implemented212\n"); }
void dsyevx_() { fprintf(stderr, "not implemented213\n"); }
void dspev_() { fprintf(stderr, "not implemented214\n"); }
void dspevx_() { fprintf(stderr, "not implemented215\n"); }
void dsbev_() { fprintf(stderr, "not implemented 216\n"); }
void dsbevx_() { fprintf(stderr, "not implemented 217\n"); }
void dsyevd_() { fprintf(stderr, "not implemented 218\n"); }
void dspevd_() { fprintf(stderr, "not implemented 219\n"); }
void dsbevd_() { fprintf(stderr, "not implemented 220\n"); }
void dsyevr_() { fprintf(stderr, "not implemented 221\n"); }
void dormbr_() { fprintf(stderr, "not implemented 222\n"); }
void dtrsyl_() { fprintf(stderr, "not implemented 223\n"); }
void dtrexc_() { fprintf(stderr, "not implemented 224\n"); }
void dtrsna_() { fprintf(stderr, "not implemented 225\n"); }
void dtrsen_() { fprintf(stderr, "not implemented 226\n"); }
void dgeevx_() { fprintf(stderr, "not implemented 227\n"); }
void dgeesx_() { fprintf(stderr, "not implemented 228\n"); }
void dggsvd_() { fprintf(stderr, "not implemented 229\n"); }
void dggsvp_() { fprintf(stderr, "not implemented 230\n"); }
void dtgsja_() { fprintf(stderr, "not implemented 231\n"); }
void dggglm_() { fprintf(stderr, "not implemented 232\n"); }
void dgglse_() { fprintf(stderr, "not implemented 233\n"); }
void dggqrf_() { fprintf(stderr, "not implemented 234\n"); }
void dggrqf_() { fprintf(stderr, "not implemented 235\n"); }
void dtgexc_() { fprintf(stderr, "not implemented 236\n"); }
void dtgsen_() { fprintf(stderr, "not implemented 237\n"); }
void dtgsna_() { fprintf(stderr, "not implemented 238\n"); }
void dtgsyl_() { fprintf(stderr, "not implemented 239\n"); }
void dormtr_() { fprintf(stderr, "not implemented 240\n"); }
void dopmtr_() { fprintf(stderr, "not implemented 241\n"); }
void dlaln2_() { fprintf(stderr, "not implemented 242\n"); }
void dlasy2_() { fprintf(stderr, "not implemented 243\n"); }
void dlanv2_() { fprintf(stderr, "not implemented 244\n"); }
void dlaexc_() { fprintf(stderr, "not implemented 245\n"); }
void dlaqtr_() { fprintf(stderr, "not implemented 246\n"); }

} // extern "C"

#endif // DUMMY

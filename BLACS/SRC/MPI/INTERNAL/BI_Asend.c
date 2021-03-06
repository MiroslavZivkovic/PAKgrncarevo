#include "Bdef.h"

void BI_Asend(BLACSCONTEXT *ctxt, int dest, int msgid, BLACBUFF *bp)
{
   int i, info, errclass;

   BI_MPI_Isend(bp->Buff, bp->N, bp->dtype, dest, msgid, ctxt->scp->comm,
                &bp->Aops[bp->nAops], info);
   while(info != BI_MPI_SUCCESS)
   {
      BI_MPI_Error_class(info, &errclass, i);
      if ( (errclass != BI_MPI_ERR_UNKNOWN) && (errclass != BI_MPI_ERR_OTHER) &&
           (errclass != BI_MPI_ERR_INTERN) )
      {
	  Mmpierror(info, "MPI_Isend", ctxt, __LINE__, __FILE__);
	  BI_BlacsErr(BI_ContxtNum(ctxt), __LINE__, __FILE__,
		      "MPI error %d on call to MPI_Isend", info);
      }
#if (BlacsDebugLvl > 0)
      else BI_BlacsWarn(BI_ContxtNum(ctxt), __LINE__, __FILE__,
"MPI error %d assumed to mean out of non-blocking resources on call to MPI_Isend",
                        info);
#endif
      BI_MPI_Isend(bp->Buff, bp->N, bp->dtype, dest, msgid, ctxt->scp->comm,
                   &bp->Aops[bp->nAops], info);
   }
   bp->nAops++;
}

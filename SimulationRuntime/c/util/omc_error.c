/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-2014, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF THE BSD NEW LICENSE OR THE
 * GPL VERSION 3 LICENSE OR THE OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the OSMC (Open Source Modelica Consortium)
 * Public License (OSMC-PL) are obtained from OSMC, either from the above
 * address, from the URLs: http://www.openmodelica.org or
 * http://www.ida.liu.se/projects/OpenModelica, and in the OpenModelica
 * distribution. GNU version 3 is obtained from:
 * http://www.gnu.org/copyleft/gpl.html. The New BSD License is obtained from:
 * http://www.opensource.org/licenses/BSD-3-Clause.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, EXCEPT AS
 * EXPRESSLY SET FORTH IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE
 * CONDITIONS OF OSMC-PL.
 *
 */

#include "setjmp.h"
#include <stdio.h>
#include "omc_error.h"
/* For MMC_THROW, so we can end this thing */
#include "meta/meta_modelica.h"

const int firstOMCErrorStream = 3;

const char *LOG_STREAM_NAME[SIM_LOG_MAX] = {
  "LOG_UNKNOWN",
  "stdout",
  "assert",

  "LOG_DASSL",
  "LOG_DASSL_STATES",
  "LOG_DEBUG",
  "LOG_DSS",
  "LOG_DSS_JAC",
  "LOG_DT",
  "LOG_EVENTS",
  "LOG_EVENTS_V",
  "LOG_INIT",
  "LOG_IPOPT",
  "LOG_IPOPT_FULL",
  "LOG_IPOPT_JAC",
  "LOG_IPOPT_HESSE",
  "LOG_IPOPT_ERROR",
  "LOG_JAC",
  "LOG_LS",
  "LOG_LS_V",
  "LOG_NLS",
  "LOG_NLS_V",
  "LOG_NLS_HOMOTOPY",
  "LOG_NLS_JAC",
  "LOG_NLS_JAC_TEST",
  "LOG_NLS_RES",
  "LOG_NLS_EXTRAPOLATE",
  "LOG_RES_INIT",
  "LOG_SIMULATION",
  "LOG_SOLVER",
  "LOG_SOLVER_CONTEXT",
  "LOG_SOTI",
  "LOG_STATS",
  "LOG_STATS_V",
#ifdef USE_DEBUG_TRACE
  "LOG_TRACE",
#endif
  "LOG_UTIL",
  "LOG_ZEROCROSSINGS",
};

const char *LOG_STREAM_DESC[SIM_LOG_MAX] = {
  "unknown",
  "this stream is always active",                       /* LOG_STDOUT */
  "this stream is always active",                       /* LOG_ASSERT */

  "additional information about dassl solver",          /* LOG_DASSL */
  "outputs the states at every dassl call",             /* LOG_DASSL_STATES */
  "additional debug information",                       /* LOG_DEBUG */
  "outputs information about dynamic state selection",  /* LOG_DSS */
  "outputs jacobian of the dynamic state selection",    /* LOG_DSS_JAC */
  "additional information about dynamic tearing",       /* LOG_DT */
  "additional information during event iteration",      /* LOG_EVENTS */
  "verbose logging of event system",                    /* LOG_EVENTS_V */
  "additional information during initialization",       /* LOG_INIT */
  "information from Ipopt",                             /* LOG_IPOPT */
  "more information from Ipopt",                        /* LOG_IPOPT_FULL*/
  "check jacobian matrix with Ipopt",                   /* LOG_IPOPT_JAC*/
  "check hessian matrix with Ipopt",                    /* LOG_IPOPT_HESSE*/
  "print max error in the optimization",                /* LOG_IPOPT_ERROR*/
  "outputs the jacobian matrix used by dassl",          /* LOG_JAC */
  "logging for linear systems",                         /* LOG_LS */
  "verbose logging of linear systems",                  /* LOG_LS_V */
  "logging for nonlinear systems",                      /* LOG_NLS */
  "verbose logging of nonlinear systems",               /* LOG_NLS_V */
  "logging of homotopy solver for nonlinear systems",   /* LOG_NLS_HOMOTOPY */
  "outputs the jacobian of nonlinear systems",          /* LOG_NLS_JAC */
  "tests the analytical jacobian of nonlinear systems", /* LOG_NLS_JAC_TEST */
  "outputs every evaluation of the residual function",  /* LOG_NLS_RES */
  "outputs debug information about extrapolate process",/* LOG_NLS_EXTRAPOLATE */
  "outputs residuals of the initialization",            /* LOG_RES_INIT */
  "additional information about simulation process",    /* LOG_SIMULATION */
  "additional information about solver process",        /* LOG_SOLVER */
  "context information during the solver process",      /* LOG_SOLVER_CONTEXT" */
  "final solution of the initialization",               /* LOG_SOTI */
  "additional statistics about timer/events/solver",    /* LOG_STATS */
  "additional statistics for LOG_STATS",                /* LOG_STATS_V */
#ifdef USE_DEBUG_TRACE
  "enables additional output to trace call stack",      /* LOG_TRACE */
#endif
  "???",                                                /* LOG_UTIL*/
  "additional information about the zerocrossings"      /* LOG_ZEROCROSSINGS */
};

static const char *LOG_TYPE_DESC[LOG_TYPE_MAX] = {
  "unknown",
  "info",
  "warning",
  "error",
  "assert",
  "debug"
};

int useStream[SIM_LOG_MAX];
int level[SIM_LOG_MAX];
int lastType[SIM_LOG_MAX];
int lastStream = LOG_UNKNOWN;
int showAllWarnings = 0;

#ifdef USE_DEBUG_TRACE
  int DEBUG_TRACE_PUSH_HELPER(const char* pFnc, const char* pFile, const long ln){if(useStream[LOG_TRACE]) printf("TRACE: push %s (%s:%d)\n", pFnc, pFile, ln); return 0;}
  int DEBUG_TRACE_POP_HELPER(int traceID){if(useStream[LOG_TRACE]) printf("TRACE: pop\n"); return 0;}
#endif

void initDumpSystem()
{
  int i;

  for(i=0; i<SIM_LOG_MAX; ++i)
  {
    useStream[i] = 0;
    level[i] = 0;
    lastType[i] = 0;
  }

  useStream[LOG_STDOUT] = 1;
  useStream[LOG_ASSERT] = 1;
}

void printInfo(FILE *stream, FILE_INFO info)
{
  fprintf(stream, "[%s:%d:%d-%d:%d:%s]", info.filename, info.lineStart, info.colStart, info.lineEnd, info.colEnd, info.readonly ? "readonly" : "writable");
}

void omc_assert_function(threadData_t* threadData, FILE_INFO info, const char *msg, ...)
{
  va_list ap;
  va_start(ap,msg);
  printInfo(stderr, info);
  fputs("Modelica Assert: ", stderr);
  vfprintf(stderr,msg,ap);
  fputs("!\n", stderr);
  va_end(ap);
  fflush(NULL);
  if (threadData) {
    MMC_THROW_INTERNAL();
  } else {
    MMC_THROW();
  }
}

void omc_assert_warning_function(FILE_INFO info, const char *msg, ...)
{
  va_list ap;
  va_start(ap,msg);
  printInfo(stderr, info);
  fputs("Warning, assertion triggered: ", stderr);
  vfprintf(stderr,msg,ap);
  fputs("!\n", stderr);
  va_end(ap);
  fflush(NULL);
}

void omc_throw_function(threadData_t *threadData)
{
  if (threadData) {
    MMC_THROW_INTERNAL();
  } else {
    MMC_THROW();
  }
}

void omc_terminate_function(FILE_INFO info, const char *msg, ...)
{
  va_list ap;
  va_start(ap,msg);
  printInfo(stderr, info);
  fputs("Modelica Terminate: ", stderr);
  vfprintf(stderr,msg,ap);
  fputs("!\n", stderr);
  va_end(ap);
  fflush(NULL);
  MMC_THROW();
}

static void printEscapedXML(const char *msg)
{
  while (*msg) {
    if (*msg == '&') fputs("&amp;", stdout);
    else if (*msg == '<') fputs("&lt;", stdout);
    else if (*msg == '>') fputs("&gt;", stdout);
    else if (*msg == '"') fputs("&quot;", stdout);
    else fputc(*msg, stdout);
    msg++;
  }
}

void messageText(int type, int stream, int indentNext, char *msg, int subline, const int *indexes)
{
  int i;
  int len;

  printf("%-17s | ", (subline || (lastStream == stream && level[stream] > 0)) ? "|" : LOG_STREAM_NAME[stream]);
  printf("%-7s | ", (subline || (lastStream == stream && lastType[stream] == type && level[stream] > 0)) ? "|" : LOG_TYPE_DESC[type]);
  lastType[stream] = type;
  lastStream = stream;

  for(i=0; i<level[stream]; ++i)
      printf("| ");

  for(i=0; msg[i]; i++)
  {
    if(msg[i] == '\n')
    {
      msg[i] = '\0';
      printf("%s\n", msg);
      if (msg[i+1]) {
        messageText(type, stream, 0, &msg[i+1], 1, indexes);
      }
      return;
    }
  }

  len = strlen(msg);
  if (len>0 && msg[len-1]=='\n') {
    printf("%s", msg);
  } else {
    printf("%s\n", msg);
  }
  fflush(NULL);
  if (indentNext) level[stream]++;
}

void messageXML(int type, int stream, int indentNext, char *msg, int subline, const int *indexes)
{
  printf("<message stream=\"%s\" type=\"%s\" text=\"", LOG_STREAM_NAME[stream], LOG_TYPE_DESC[type]);
  printEscapedXML(msg);
  if (indexes) {
    int i;
    printf("\">\n");
    for (i=1; i<=*indexes; i++) {
      printf("<used index=\"%d\" />\n", indexes[i]);
    }
    if (!indentNext) {
      fputs("</message>\n",stdout);
    }
  } else {
    fputs(indentNext ? "\">\n" : "\" />\n", stdout);
  }
  fflush(stdout);
}

static void messageCloseText(int stream)
{
  if(ACTIVE_STREAM(stream))
    level[stream]--;
}

static void messageCloseXML(int stream)
{
  if(ACTIVE_STREAM(stream))
  {
    fputs("</message>\n", stdout);
    fflush(stdout);
  }
}

static void messageCloseTextWarning(int stream)
{
  if(ACTIVE_WARNING_STREAM(stream))
    level[stream]--;
}

static void messageCloseXMLWarning(int stream)
{
  if(ACTIVE_WARNING_STREAM(stream))
  {
    fputs("</message>\n", stdout);
    fflush(stdout);
  }
}

static void (*messageFunction)(int type, int stream, int indentNext, char *msg, int subline, const int *indexes) = messageText;
void (*messageClose)(int stream) = messageCloseText;
void (*messageCloseWarning)(int stream) = messageCloseTextWarning;

void setStreamPrintXML(int isXML)
{
  if (isXML) {
    messageFunction = messageXML;
    messageClose = messageCloseXML;
    messageCloseWarning = messageCloseXMLWarning;
  } else {
    messageFunction = messageText;
    messageClose = messageCloseText;
    messageCloseWarning = messageCloseTextWarning;
  }
}

#define SIZE_LOG_BUFFER 2048
void va_infoStreamPrint(int stream, int indentNext, const char *format, va_list args)
{
  if (useStream[stream]) {
    char logBuffer[SIZE_LOG_BUFFER];
    vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
    messageFunction(LOG_TYPE_INFO, stream, indentNext, logBuffer, 0, NULL);
  }
}

void infoStreamPrintWithEquationIndexes(int stream, int indentNext, const int *indexes, const char *format, ...)
{
  if (useStream[stream]) {
    char logBuffer[SIZE_LOG_BUFFER];
    va_list args;
    va_start(args, format);
    vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
    va_end(args);
    messageFunction(LOG_TYPE_INFO, stream, indentNext, logBuffer, 0, indexes);
  }
}

void infoStreamPrint(int stream, int indentNext, const char *format, ...)
{
  if (useStream[stream]) {
    char logBuffer[SIZE_LOG_BUFFER];
    va_list args;
    va_start(args, format);
    vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
    va_end(args);
    messageFunction(LOG_TYPE_INFO, stream, indentNext, logBuffer, 0, NULL);
  }
}

void warningStreamPrintWithEquationIndexes(int stream, int indentNext, const int *indexes, const char *format, ...)
{
  if (ACTIVE_WARNING_STREAM(stream)) {
    char logBuffer[SIZE_LOG_BUFFER];
    va_list args;
    va_start(args, format);
    vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
    va_end(args);
    messageFunction(LOG_TYPE_WARNING, stream, indentNext, logBuffer, 0, indexes);
  }
}

void warningStreamPrint(int stream, int indentNext, const char *format, ...)
{
  if (ACTIVE_WARNING_STREAM(stream)) {
    char logBuffer[SIZE_LOG_BUFFER];
    va_list args;
    va_start(args, format);
    vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
    va_end(args);
    messageFunction(LOG_TYPE_WARNING, stream, indentNext, logBuffer, 0, NULL);
  }
}

void va_warningStreamPrintWithEquationIndexes(int stream, int indentNext, const int *indexes, const char *format, va_list args)
{
  if (ACTIVE_WARNING_STREAM(stream)) {
    char logBuffer[SIZE_LOG_BUFFER];
    vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
    messageFunction(LOG_TYPE_WARNING, stream, indentNext, logBuffer, 0, indexes);
  }
}

void va_warningStreamPrint(int stream, int indentNext, const char *format, va_list args)
{
  if (ACTIVE_WARNING_STREAM(stream)) {
    char logBuffer[SIZE_LOG_BUFFER];
    vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
    messageFunction(LOG_TYPE_WARNING, stream, indentNext, logBuffer, 0, NULL);
  }
}

void errorStreamPrint(int stream, int indentNext, const char *format, ...)
{
  char logBuffer[SIZE_LOG_BUFFER];
  va_list args;
  va_start(args, format);
  vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
  va_end(args);
  messageFunction(LOG_TYPE_ERROR, stream, indentNext, logBuffer, 0, NULL);
}

void va_errorStreamPrint(int stream, int indentNext, const char *format, va_list args)
{
  char logBuffer[SIZE_LOG_BUFFER];
  vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
  messageFunction(LOG_TYPE_ERROR, stream, indentNext, logBuffer, 0, NULL);
}

void va_errorStreamPrintWithEquationIndexes(int stream, int indentNext, const int *indexes, const char *format, va_list args)
{

  char logBuffer[SIZE_LOG_BUFFER];
  vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
  messageFunction(LOG_TYPE_ERROR, stream, indentNext, logBuffer, 0, indexes);
}

#ifdef USE_DEBUG_OUTPUT
void debugStreamPrint(int stream, int indentNext, const char *format, ...)
{
  if (useStream[stream]) {
    char logBuffer[SIZE_LOG_BUFFER];
    va_list args;
    va_start(args, format);
    vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
    va_end(args);
    messageFunction(LOG_TYPE_DEBUG, stream, indentNext, logBuffer, 0, NULL);
  }
}

void debugStreamPrintWithEquationIndexes(int stream, int indentNext, const int *indexes, const char *format, ...)
{
  if (useStream[stream]) {
    char logBuffer[SIZE_LOG_BUFFER];
    va_list args;
    va_start(args, format);
    vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
    va_end(args);
    messageFunction(LOG_TYPE_DEBUG, stream, indentNext, logBuffer, 0, indexes);
  }
}
#endif

static inline jmp_buf* getBestJumpBuffer(threadData_t *threadData)
{
  switch (threadData->currentErrorStage) {
  case ERROR_EVENTSEARCH:
  case ERROR_SIMULATION:
  case ERROR_NONLINEARSOLVER:
  case ERROR_INTEGRATOR:
  case ERROR_OPTIMIZE:
#ifndef OMC_EMCC
    return threadData->simulationJumpBuffer;
#endif
  case ERROR_EVENTHANDLING:
  default:
    return threadData->globalJumpBuffer ? threadData->globalJumpBuffer : threadData->mmc_jumper;
  }
}

void va_throwStreamPrint(threadData_t *threadData, const char *format, va_list args)
{
  char logBuffer[SIZE_LOG_BUFFER];
  vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
  messageFunction(LOG_TYPE_DEBUG, LOG_ASSERT, 0, logBuffer, 0, NULL);
  threadData = threadData ? threadData : (threadData_t*)pthread_getspecific(mmc_thread_data_key);
  longjmp(*getBestJumpBuffer(threadData), 1);
}

void throwStreamPrint(threadData_t *threadData, const char *format, ...)
{
  va_list args;
  va_start(args, format);
  va_throwStreamPrint(threadData, format, args);
  va_end(args);
}

void throwStreamPrintWithEquationIndexes(threadData_t *threadData, const int *indexes, const char *format, ...)
{
  char logBuffer[SIZE_LOG_BUFFER];
  va_list args;
  va_start(args, format);
  vsnprintf(logBuffer, SIZE_LOG_BUFFER, format, args);
  va_end(args);
  messageFunction(LOG_TYPE_DEBUG, LOG_ASSERT, 0, logBuffer, 0, indexes);
  threadData = threadData ? threadData : (threadData_t*)pthread_getspecific(mmc_thread_data_key);
  longjmp(*getBestJumpBuffer(threadData), 1);
}

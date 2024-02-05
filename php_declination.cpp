#define WIN32_LEAN_AND_MEAN

#include <windows.h>

#include "zend_config.w32.h"
#include "php.h"
#include "ext/standard/info.h"

#include "GeoMagLib.h"

extern zend_module_entry declination_module_entry;
#define phpext_declination_ptr &declination_module_entry

#define PHP_TEST_VERSION "1.0"

#if defined(ZTS) && defined(COMPILE_DL_TEST)
ZEND_TSRMLS_CACHE_EXTERN()
#endif

ZEND_BEGIN_MODULE_GLOBALS(declination)
zend_long someGlobal;
ZEND_END_MODULE_GLOBALS(declination)

ZEND_EXTERN_MODULE_GLOBALS(declination)

#define TEST_G(v) ZEND_MODULE_GLOBALS_ACCESSOR(declination, v)

CGeoMagnetic* geoMagnetic = nullptr;

double sdate = 2021.15;

char projection = 'M';

double alt = 0.0;

double longitude = -81.440772;
double latitude = 39.342776;

char unit = 'F';

#pragma comment(lib, "php8ts.lib")

#pragma comment(lib, "GeoMagLib.lib")

/* For compatibility with older PHP versions */
#ifndef ZEND_PARSE_PARAMETERS_NONE
#define ZEND_PARSE_PARAMETERS_NONE() \
	ZEND_PARSE_PARAMETERS_START(0, 0) \
	ZEND_PARSE_PARAMETERS_END()
#endif

ZEND_DECLARE_MODULE_GLOBALS(declination)

PHP_INI_BEGIN()
	STD_PHP_INI_ENTRY("declination.someGlobal", "1", PHP_INI_ALL, OnUpdateLong, someGlobal,	zend_declination_globals, declination_globals)
PHP_INI_END()

PHP_FUNCTION(getDeclination)
{
	ZEND_PARSE_PARAMETERS_START(3, 3)
		Z_PARAM_DOUBLE(latitude)
		Z_PARAM_DOUBLE(longitude)
		Z_PARAM_DOUBLE(sdate)
	ZEND_PARSE_PARAMETERS_END();

	geoMagnetic = new CGeoMagnetic(sdate, projection);

	geoMagnetic->CalculateFieldElements(latitude, longitude, alt, unit);

	double d = geoMagnetic->GeoMagneticElements->Decl;

	delete geoMagnetic;

	RETURN_DOUBLE(d);
}

static PHP_GINIT_FUNCTION(declination)
{
#if defined(COMPILE_DL_BCMATH) && defined(ZTS)
	ZEND_TSRMLS_CACHE_UPDATE();
#endif
	declination_globals->someGlobal = 1;
}

PHP_MINIT_FUNCTION(declination)
{
	REGISTER_INI_ENTRIES();

	return SUCCESS;
}

PHP_MINFO_FUNCTION(declination)
{
	php_info_print_table_start();
	php_info_print_table_row(2, "Magnetic Declination Support", "enabled");
	php_info_print_table_row(2, "The World Magnetic Model", "https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml");
	php_info_print_table_row(2, "Function", "double decl = getDeclination(latitude, longitude, yyyy.ddd);");
	php_info_print_table_row(2, "Version", "1.0");
	php_info_print_table_end();
}

ZEND_BEGIN_ARG_INFO(arginfo_getDeclination, 0)
ZEND_ARG_INFO(0, latitude)
ZEND_ARG_INFO(0, longitude)
ZEND_ARG_INFO(0, sdate)
ZEND_END_ARG_INFO()

static const zend_function_entry declination_functions[] =
{
	PHP_FE(getDeclination,	arginfo_getDeclination)
	PHP_FE_END
};

zend_module_entry declination_module_entry =
{
	STANDARD_MODULE_HEADER,
	"declination",						/* Extension name */
	declination_functions,				/* zend_function_entry */
	PHP_MINIT(declination),				/* PHP_MINIT - Module initialization */
	NULL,								/* PHP_MSHUTDOWN - Module shutdown */
	NULL,								/* PHP_RINIT - Request initialization */
	NULL,								/* PHP_RSHUTDOWN - Request shutdown */
	PHP_MINFO(declination),				/* PHP_MINFO - Module info */
	PHP_TEST_VERSION,					/* Version */
	PHP_MODULE_GLOBALS(declination),	/* Module globals */
	PHP_GINIT(declination),				/* PHP_GINIT - Globals initialization */
	NULL,								/* PHP_GSHUTDOWN - Globals shutdown */
	NULL,
	STANDARD_MODULE_PROPERTIES_EX
};

#ifdef COMPILE_DL_TEST
# ifdef ZTS
ZEND_TSRMLS_CACHE_DEFINE()
# endif
ZEND_GET_MODULE(declination)
#endif
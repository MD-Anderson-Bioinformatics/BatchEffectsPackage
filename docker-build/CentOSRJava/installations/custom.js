//==============================================//
// Standard link out file for NG-CHMs           //
//==============================================//
linkouts.setVersion('2020-04-07');

// Scatter Plot plugin:
linkouts.addPanePlugin ({
	name: 'ScatterPlot',
	helpText: 'Creates a two-dimensional scatter plot',
	params: {},
	src: 'https://www.ngchm.net/Plugins/ScatterPlot/index.html'
});

//==============================================//
// Biology base plugin                          //
//==============================================//

// PathwayMapper modified to be an NG-CHM plugin:
linkouts.addPanePlugin ({
	name: 'PathwayMapper',
	helpText: 'Pathway Mapper',
	params: {},
	src: 'https://www.ngchm.net/Plugins/PathwayMapper/index.html'
});

(function(linkouts) {

	/*BEGIN SAMPLE Linkouts to the Hamburger Menu
	//TODO: replace with "real" hamburger linkout(s)
	//Using default icon
	linkouts.addHamburgerLinkout({name: "menuLink2", label: "Execute Menu Link 1", action: "linkoutHelp();"});
	//Using different icon from server
	linkouts.addHamburgerLinkout({name: "menuLink3", label: "Execute Menu Link 2", icon: "images/redX.png", action: "linkoutHelp();"}); 
	//Using external icon from the web
	linkouts.addHamburgerLinkout({name: "menuLink4", label: "Execute Menu Link 3", icon: "http://amigo.geneontology.org/static/images/go-logo-icon.small.png", action: "openWidgetHelp();"});
	//No label provided
	linkouts.addHamburgerLinkout({name: "menuLink5", action: "openWidgetHelp();"});
	//No action provided
	linkouts.addHamburgerLinkout({name: "menuLink6"});
	//END SAMPLE Linkouts to the Hamburger Menu*/
 
	//Add Heat Map Linkouts
	linkouts.addSubtype ("bio.gene.hugo", "bio.pubmed");
    linkouts.addSubtype ("bio.mirna", "bio.pubmed");
    linkouts.addSubtype ("bio.gene.entrez", "bio.gene.entrezid");

    linkouts.describeTypes([
        { typeName: "bio.gene.hugo",
          displayName: "HUGO gene symbol",
          description: "The official HUGO symbol for a gene.",
          examples: "TP53",
          format: "Latin letters and Arabic numerals, usually without punctuation."
        },
        { typeName: "bio.go",
          displayName: "Gene Ontology accession identifier",
          description: "A Gene Ontology (GO) accession identifier.",
          examples: "GO:0007015",
          format: "Letters G and O, a colon, and a seven-digit number."
        },
        { typeName: "bio.transcriptid",
          displayName: "Entrez transcription identifier",
          description: "An Entrez transcription identifier",
          format: "ENST###########"
        },
        { typeName: "bio.mirna",
          displayName: "miRNA identifier",
          description: "An miRNA identifier"
        },
        { typeName: "bio.pathway.msigdb.name",
          displayName: "msig DB pathway name",
          description: "A Pathway Name as defined by the msig DB"
        },
        { typeName: "bio.gene.entrezid",
          displayName: "NCBI Entrez Gene Identifier",
          description: "Numeric Id representing a gene as defined by NCBI"
        },
        { typeName: "bio.protein.uniprotid",
          displayName: "Uniprot protein identifier",
          description: "A Uniprot protein identifier."
        },
        { typeName: "bio.geo.acc",
          displayName: "GEO accession id",
          description: "A Gene Expression Omnibus (GEO) accession id"
        },
        { typeName: "bio.pubmed",
          displayName: "PubMed identifier",
          description: "A PubMed Identifier, (Integer ID)",
          format: "an integer"
        },
        { typeName: "bio.transcript.ensembl",
          displayName: "Ensembl transcript identifier",
          description: "An Ensembl transcript identifier"
        },
        { typeName: "scholar",
          displayName: "Google Scholar search term",
          description: "A text string used to search for something in Google Scholar"
        },
        { typeName: "search",
          displayName: "Generic search term",
          description: "A text string used to search for something"
        }
    ]);
}) (linkouts);

/*BEGIN Sample Hamburger Linkout Functions
 TODO: replace with "real" hamburger linkout function(s)
function openWidgetHelp () {
	NgChm.UHM.widgetHelp();
}

function linkoutHelp () {
	NgChm.UHM.openLinkoutHelp();
}
//END Sample Hamburger Linkout Functions*/
 
//==============================================//
// Amigo plugin                                 //
//==============================================//
(function(linkouts) {

    function openAmigo (names) {
	var goid = names[0];
	linkouts.openUrl("http://amigo.geneontology.org/amigo/term/" + goid, "GeneOntology");
    };

    linkouts.addPlugin({
        name: "Amigo",
	description: "Adds linkouts to Amigo gene ontology database.",
	version: "0.1.0",
	site: "http://amigo.geneontology.org/amigo",
	logo: "http://amigo.geneontology.org/static/images/go-logo-icon.small.png",
	linkouts: [
	    { menuEntry: "View Amigo", typeName: "bio.go", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openAmigo }
	]
    });
}) (linkouts);

//==============================================//
// BioGPS plugin                                //
//==============================================//
(function(linkouts) {

    function searchBioGPS (names) {
	var gname = names[0];
	linkouts.openUrl("http://biogps.org/search/?q=" + gname, "bioGPS", { noframe: true });
    }

    linkouts.addPlugin({
        name: "BioGPS",
        description: "Adds linkouts to the BioGPS gene annotation portal.",
        version: "0.1.0",
        site: "http://biogps.org/",
        logo: "http://biogps.org/assets/img2/biogps_logo.png",
        linkouts: [
	    { menuEntry: "View bioGPS", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: searchBioGPS }
        ]
    });
}) (linkouts);

//==============================================//
// Cancer Digital Slide Archive plugin          //
//==============================================//
(function(linkouts) {
    function openSlideArchive (ids) {
        var part = ids[0].substr(0,12);
        linkouts.openUrl("https://cancer.digitalslidearchive.org/#!/CDSA/byPatientID/" + part, "SlideArchive");
    }
    linkouts.addPlugin({
        name: "Cancer Digital Slide Archive",
        site: "https://cancer.digitalslidearchive.org/",
        logo: "https://cancer.digitalslidearchive.org/img/CDSA_Slide_50.png",
        description: "Adds linkouts to the Cancer Digitial Slide Archive of TCGA digital slide images.",
        version: "0.3.0",
        linkouts: [
            { menuEntry: "View SlideArchive", typeName: "bio.tcga.barcode.sample", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openSlideArchive }
        ]});
})(linkouts);

//==============================================//
// cBioPortal plugin                            //
//==============================================//
(function (linkouts) {

    linkouts.describeTypes([
        { typeName: "bio.cbioportal.sampleid",
          displayName: "cBioPortal Sample Identifier",
          description: "An identifier that refers to a specific sample in cBioPortal."
        }
    ]);

    function openCBIOGenes ( genes) {
	var studyid = linkouts.getAttribute("cbioportal.study");
	linkouts.openUrl("http://www.cbioportal.org/ln?cancer_study_id=" + studyid + "&gene_list=" + genes.join("+"), "cbio");
    }

    function openCBIOSamp (ids) {
        var studyid = linkouts.getAttribute("cbioportal.study");
        linkouts.openUrl("http://www.cbioportal.org/case.do?cancer_study_id=" + studyid + "&case_id=" + ids[0], "cbio");
    }

    function openCBIOSampTCGA (ids) {
        var part = ids[0].substring(0,12);
        openCBIOSamp ([part]);
    }

    linkouts.addPlugin({
        name: "cBioPortal",
        description: "Adds linkouts to the cBioPortal for Cancer Genomics.",
        version: "0.1.0",
        logo: "http://www.cbioportal.org/images/cbioportal_logo.png",
        site: "http://www.cbioportal.org/",
        attributes: [
            { name: "cbioportal.study", description: "cBioPortal study identifier" }
        ],
        linkouts: [
            { menuEntry: "View in cBio Portal", typeName: "bio.cbioportal.sampleid", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openCBIOSamp, attributes: ["cbioportal.study"] },
            { menuEntry: "View in cBio Portal", typeName: "bio.tcga.barcode.sample", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openCBIOSampTCGA, attributes: ["cbioportal.study"] },
            { menuEntry: "View in cBio Portal", typeName: "bio.gene.hugo", selectMode: linkouts.MULTI_SELECT, linkoutFn: openCBIOGenes, attributes: ["cbioportal.study"] }
        ]});
})(linkouts);

//==============================================//
// CivicDB plugin                               //
//==============================================//
(function(linkouts) {
    function openCivicDB (ids) {
	linkouts.openUrl("https://civicdb.org/links/entrez_id/" + ids[0]);
    }
    linkouts.addPlugin({
        name: "CIViC",
        description: "Adds linkouts to the CIViC Mutation Database",
        version: "0.1.0",
        site: "https://civic.genome.wustl.edu/home",
        logo: "https://civic.genome.wustl.edu/assets/images/CIViC_logo.png",
        linkouts: [
            { menuEntry: "View CIViC", typeName: "bio.gene.entrezid", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openCivicDB }
	]
    });
}) (linkouts);

//==============================================//
// COSMIC plugin                                //
//==============================================//
(function(linkouts) {

    function openCosmicGene (names) {
	var gname = names[0];
	linkouts.openUrl("http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=" + gname, "Cosmic");
    }

    linkouts.addPlugin({
        name: "COSMIC",
        description: "Adds linkouts to the Catalogue of somatic mutations in cancer (COSMIC).",
        version: "0.1.0",
        site: "http://cancer.sanger.ac.uk/cosmic",
        logo: "http://cancer.sanger.ac.uk/cancergenome/gfx/logo_cosmic.png",
        linkouts: [
	    { menuEntry: "View Cosmic", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openCosmicGene }
        ]
    });
}) (linkouts);

//==============================================//
// Decipher plugin                              //
//==============================================//
(function(linkouts) {

    function openDecipher (names) {
	var gname = names[0];
	linkouts.openUrl("https://decipher.sanger.ac.uk/search?q=" + gname, "Decipher", { noframe: true });
    }

    linkouts.addPlugin({
        name: "Decipher",
        description: "Adds linkouts to the Decipher database",
        version: "0.1.0",
        site: "https://decipher.sanger.ac.uk/",
        logo: "https://decipher.sanger.ac.uk/img/decipher-logo.png",
        linkouts: [
            { menuEntry: "View Decipher", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openDecipher }
        ]
    });
}) (linkouts);

//==============================================//
// DepMap plugin                                //
//==============================================//
(function(linkouts) {
    function openDepMap (ids) {
	linkouts.openUrl("https://depmap.org/portal/gene/" + ids[0] + "?tab=overview");
    }
    linkouts.addPlugin({
        name: "DepMap",
        description: "Adds linkouts to the DepMap Portal",
        version: "0.1.0",
        site: "https://depmap.org/portal/",
        logo: "https://depmap.org/portal/static/img/nav_footer/banner_logo_depmapportal.svg",
        linkouts: [
            { menuEntry: "View depmap", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openDepMap }
	]
    });
}) (linkouts);

//==============================================//
// Ensembl plugin                               //
//==============================================//
(function(linkouts) {

    function searchEnsemblForGene (names) {
	const gname = names[0];
	const species = linkouts.getAttribute("bio.species") || "Homo_sapiens";
	linkouts.openUrl("http://ensembl.org/Multi/psychic?site=ensembl&species="+species+"&q=" + gname, "Ensembl");
    }

    function searchEnsemblForTranscript (names) {
	const tname = names[0];
	const species = linkouts.getAttribute("bio.species") || "Homo_sapiens";
	linkouts.openUrl("http://ensembl.org/Multi/psychic?site=ensembl&species="+species+"&q=" + tname, "Ensembl");
    }

    linkouts.addPlugin({
        name: "Ensembl",
	description: "Adds linkouts to Ensembl genome browser.",
	version: "0.1.0",
	site: "http://www.ensembl.org/index.html",
	logo: "http://www.sanger.ac.uk/sites/default/files/ensembl_1.gif",
	linkouts: [
	    // linkouts for gene symbols
            { menuEntry: "Search Ensembl", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: searchEnsemblForGene },
	    // linkouts for transcript ids
	    { menuEntry: "Search Ensembl", typeName: "bio.transcript.ensembl", selectMode: linkouts.SINGLE_SELECT, linkoutFn: searchEnsemblForTranscript }
	]
    });
}) (linkouts);

//==============================================//
// FireBrowse plugin                            //
//==============================================//
(function(linkouts) {

    function openFireBrowseGene (names) {
	var gname = names[0];
	linkouts.openUrl("http://firebrowse.org/viewGene.html?gene=" + gname + "&search=" + gname, "FireBrowse");
    }

    linkouts.addPlugin({
        name: "FireBrowse",
	description: "Adds linkouts to FireBrowse.",
	version: "0.1.0",
	site: "http://firebrowse.org/",
	logo: "http://firebrowse.org/img/firebrowse.png",
	linkouts: [
            { menuEntry: "View FireBrowse", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openFireBrowseGene }
	]
    });
}) (linkouts);

//==============================================//
//GeneCards plugin                                //
//==============================================//
(function(linkouts) {

  function searchGeneCards(labels){
	    linkouts.openUrl("http://www.genecards.org/Search/Keyword?queryString=" + labels[0], "GeneCards");
  }

linkouts.addPlugin({
   name: "GeneCards",
   description: "Adds linkouts to the GeneCards Human Gene Card Database.",
   version: "0.1.2",
   site: "http://www.genecards.org/",
   logo: "http://www.genecards.org/Images/logo_genecards.png",
   linkouts: [
	    { menuEntry: "View GeneCards", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: searchGeneCards }
   ]
});
}) (linkouts);

//==============================================//
//GeneVisible plugin                                 //
//==============================================//
(function(linkouts) {

  function openGenevisiblePeptide (names) {
  	var gname = names[0];
  	linkouts.openUrl("https://genevisible.com/tissues/HS/UniProt/" + gname, "GeneVisible");
      }

      function openGenevisibleHugo (names) {
  	var gname = names[0];
  	linkouts.openUrl("https://genevisible.com/tissues/HS/Gene%20Symbol/" + gname, "GeneVisible");
      }


linkouts.addPlugin({
   name: "GeneVisible",
	description: "Adds linkouts to the GeneVisible portal.",
	version: "0.1.2",
	site: "https://genevisible.com",
	logo: "https://genevisible.com/img/genevisible_logo.svgz",
	linkouts: [
	            { menuEntry: "View GeneVisible", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openGenevisibleHugo },
	            { menuEntry: "View GeneVisible", typeName: "bio.protein.uniprotid", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openGenevisiblePeptide }
	]
});
}) (linkouts);

//==============================================//
// Google plugin                                //
//==============================================//
(function(linkouts) {

    linkouts.addSubtype ("bio.gene.hugo", "scholar");
    linkouts.addSubtype ("bio.mirna", "scholar");

    linkouts.addPlugin({
        name: "Google Scholar",
	description: "Add linkout to search Google Scholar.",
	version: "0.1.2",
	site: "https://scholar.google.com/schhp?hl=en",
	logo: "https://scholar.google.com/intl/en/scholar/images/1x/googlelogo_color_270x104dp.png",
	linkouts: [
	    { menuEntry: "Search Google Scholar", typeName: "scholar", selectMode: linkouts.MULTI_SELECT, linkoutFn: searchGoogleScholar },
	    { menuEntry: "Search Google", typeName: "search", selectMode: linkouts.MULTI_SELECT, linkoutFn: searchGoogle }
	]
    });

    function searchGoogle(selection){
	linkouts.openUrl("https://www.google.com/#q=" + selection.join("+"), { noframe: true });
    }

    function searchGoogleScholar (labels) {
        linkouts.openUrl("http://scholar.google.com/scholar?q=" + labels.join("+OR+"), "pubmed", { noframe: true });
    };

}) (linkouts);

//==============================================//
// GTEx plugin                                  //
//==============================================//
(function(linkouts) {
    function openGTExPortal (ids) {
	linkouts.openUrl("https://gtexportal.org/home/gene/" + ids[0]);
    }
    linkouts.addPlugin({
        name: "GTEx Portal",
        description: "Adds linkouts to the GTEx Portal",
        version: "0.1.0",
        site: "https://gtexportal.org/home/",
        logo: "https://gtexportal.org/img/gtex2.1a2a339c.png",
        linkouts: [
            { menuEntry: "View GTEx", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openGTExPortal }
	]
    });
}) (linkouts);

//==============================================//
//HGNC plugin                                 //
//==============================================//
(function(linkouts) {

	const working = false;
  function openHGNCGene (names) {
      var gname = names[0];
      linkouts.openUrl("http://www.genenames.org/cgi-bin/gene_symbol_report?q=data/hgnc_data.php&match=" + gname, "hgnc");
  }

if (working) {
linkouts.addPlugin({
   name: "HGNC",
	description: "Adds linkouts to HGNC portal.",
	version: "0.1.2",
	site: "http://www.genenames.org",
	logo: "https://www.genenames.org/sites/genenames.org/themes/custom/genenames/genenames_logo.png",
	linkouts: [
	    { menuEntry: "View HGNC", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openHGNCGene }
	]
});
}
}) (linkouts);

//==============================================//
// Ideogram viewer plugin                       //
//==============================================//
(function(linkouts) {

    function viewIdeogramGene(genes){
	linkouts.openUrl("http://bioinformatics.mdanderson.org/ideogramviewer/Ideogram.html?genelist1=" + genes.join(","), "Ideogram");
    }

    function viewIdeogramGene2(labels){
        var genes1 = labels.Row;
        var genes2 = labels.Column;
	linkouts.openUrl("http://bioinformatics.mdanderson.org/ideogramviewer/Ideogram.html?genelist1=" + genes1.join(",") + "&genelist2=" + genes2.join(","),
	"ideogram");
    }

    function viewIdeogramMIRNA(mir){
	mir = mir.map(function (m) { return m.substring(0, m.lastIndexOf(".MIMAT")) });
	linkouts.openUrl("http://bioinformatics.mdanderson.org/ideogramviewer/Ideogram.html?mirlist1=" + mir.join(","),
	"ideogram");
    }

    function viewIdeogramMIRNA2(labels){
	var mirs1 = labels.Row.map(function (m) { return m.substring(0, m.lastIndexOf(".MIMAT")) });
	var mirs2 = labels.Column.map(function (m) { return m.substring(0, m.lastIndexOf(".MIMAT")) });
	linkouts.openUrl("http://bioinformatics.mdanderson.org/ideogramviewer/Ideogram.html?mirlist1=" + mirs1.join(",") + "&mirlist2=" + mirs2.join(","),
	"ideogram");
    }

    linkouts.addPlugin({
        name: "Ideogram Viewer",
	description: "Adds linkouts for viewing a set of genes and/or mirs on an interactive ideogram.",
	version: "0.1.0",
	site: "http://bioinformatics.mdanderson.org/main/IdeogramViewer:Overview",
        logo: "http://bioinformatics.mdanderson.org/people/Bradley.Broom/IdeogramViewerLogo.png",
	linkouts: [
	    { menuEntry: "View ideogram", typeName: "bio.gene.hugo", selectMode: linkouts.MULTI_SELECT, linkoutFn: viewIdeogramGene },
            { menuEntry: "View ideogram", typeName: "bio.mirna", selectMode: linkouts.MULTI_SELECT, linkoutFn: viewIdeogramMIRNA }
	],
        matrixLinkouts: [
            { menuEntry: "View ideogram", typeName1: ["bio.gene.hugo"], typeName2: ["bio.gene.hugo"], selectMode: linkouts.MULTI_SELECT, linkoutFn: viewIdeogramGene2 },
            { menuEntry: "View ideogram", typeName1: ["bio.mirna"], typeName2: ["bio.mirna"], selectMode: linkouts.MULTI_SELECT, linkoutFn: viewIdeogramMIRNA2 }
        ]
    });
}) (linkouts);

//==============================================//
// MaveDB Variant Effects Database plugin       //
//==============================================//
(function(linkouts) {

    function openMavedbGene (names) {
	const gname = names[0];
	const species = linkouts.getAttribute("bio.species") || "Homo_sapiens";
	linkouts.openUrl("https://mavedb.org/search/?organism=" + species + "&search=" + gname, "MaveDB", { noframe: true });
    };

    linkouts.addPlugin({
        name: "MaveDB Variant Effects Database",
	description: "Adds linkout to search MaveDB.",
	version: "0.1.0",
	site: "https://www.mavedb.org",
	logo: "https://www.mavedb.org/static/core/mavedb/dna.png",
	linkouts: [
	    { menuEntry: "Search MaveDB", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openMavedbGene }
	]
    });
}) (linkouts);

//==============================================//
// MD Anderson Pathway Database plugin          //
//==============================================//
(function(linkouts) {

    linkouts.describeTypes([
        { typeName: "bio.mdacc.pathwayid",
          displayName: "PathwaysWeb pathway identifier",
          description: "A Pathway Id as defined by the MD Anderson PathwaysWeb System"
        }
    ]);

    function openMDACCPathwayID (names) {
	var gname = names[0];
	linkouts.openUrl("https://bioinformatics.mdanderson.org/PathwaysBrowser/pathway/latest/mdaPathwayId/" + gname, "Pathways");
    };

    function openPathwaysBrowserGO (names) {
	var goid = names[0];
	linkouts.openUrl("https://bioinformatics.mdanderson.org/PathwaysBrowser/goTerm/latest/goId/" + goid, "Geneontology");
    };

    linkouts.addPlugin({
        name: "Pathways Browser",
	description: "Adds linkouts to the MD Anderson Pathways Browser.",
	version: "0.1.0",
	site: "https://bioinformatics.mdanderson.org/PathwaysBrowser/",
	logo: "https://bioinformatics.mdanderson.org//PathwaysBrowser/mda_logo.png",
	linkouts: [
	    { menuEntry: "View PathwaysBrowser", typeName: "bio.go", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openPathwaysBrowserGO },
	    // linkouts for pathways
	    { menuEntry: "View Pathway", typeName: "bio.mdacc.pathwayid", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openMDACCPathwayID }
	]
    });
}) (linkouts);

//==============================================//
// miRBase plugin                               //
//==============================================//
(function(linkouts) {

    function viewMiRBasePage (eids) {
	var gid = eids[0];
	gid = gid.substr( gid.lastIndexOf(".MIMAT") + 1 );
	linkouts.openUrl("http://www.mirbase.org/cgi-bin/query.pl?terms=" + gid, "miRBase");
    }

    linkouts.addPlugin({
        name: "miRBase",
	description: "Adds links to miRBase.",
	version: "0.1.0",
	site: "http://www.mirbase.org/",
	logo: "http://www.mirbase.org/images/mirbase-logo-blue-web.png",
	linkouts: [
	    { menuEntry: "View miRBase Page", typeName: "bio.mirna", selectMode: linkouts.SINGLE_SELECT, linkoutFn: viewMiRBasePage }
	]
    });
}) (linkouts);

//==============================================//
// MSigDB plugin                                //
//==============================================//
(function(linkouts) {
    function openMSigDB (names) {
	var pwname = names[0];
	linkouts.openUrl("http://software.broadinstitute.org/gsea/msigdb/cards/" + pwname + ".html", "Uniprot");
    }

    linkouts.addPlugin({
        name: "MSigDB",
        description: "Adds linkouts to MSigDB pathways database",
	version: "0.1.0",
	site: "http://software.broadinstitute.org/gsea/index.jsp",
	logo: "http://software.broadinstitute.org/gsea/images/GSEA-logo.gif",
	linkouts: [
	    { menuEntry: "View MSigDB", typeName: "bio.pathway.msigdb.name", selectMode: linkouts.SINGLE_SELECT,
	      linkoutFn: openMSigDB }
	]
    });
}) (linkouts);

//==============================================//
// MuPIT plugin                                 //
//==============================================//
(function(linkouts) {

    function viewMupitG (genes) {
	genes = genes.sort().filter(function(el,i,a){return i==a.indexOf(el);});
	var glist = encodeURIComponent(genes[0]);
	linkouts.openUrl("http://mupit.icm.jhu.edu/MuPIT_Interactive?gene=" + glist, "MuPIT");
    };

    linkouts.addPlugin({
        name: "MuPIT",
        description: "Adds linkouts to MuPIT Interactive.",
        version: "0.1.0",
	site: "http://mupit.icm.jhu.edu/MuPIT_Interactive/",
	logo: "http://mupit.icm.jhu.edu/MuPIT_Interactive/images/muPITlog.gif",
	linkouts: [
            { menuEntry: "View MuPIT", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: viewMupitG }
	]
    });
}) (linkouts);

//==============================================//
// NCBI plugin                                  //
//==============================================//
(function(linkouts) {

    function openClinVar (names) {
	var gname = names[0];
	linkouts.openUrl("https://www.ncbi.nlm.nih.gov/clinvar/?term=" + gname + "%5Bgene%5D", "ClinVar");
    }

    function openNCBIGenePage (names) {
	const gname = names[0];
	const species = linkouts.getAttribute("bio.species") || "Homo_sapiens";
	linkouts.openUrl("http://www.ncbi.nlm.nih.gov/gene?term=(" + species + "%5BOrganism%5D)%20AND%20" + gname + "%5BGene%20Name%5D", "NCBI");
    }
    function openNCBIEntrezIDPage (eids) {
	var gid = eids[0];
	linkouts.openUrl("http://www.ncbi.nlm.nih.gov/gene/" + gid, "NCBI");
    }
    function searchNCBIDatabases (names) {
	const gname = names[0];
	const species = linkouts.getAttribute("bio.species") || "Homo_sapiens";
	linkouts.openUrl("http://www.ncbi.nlm.nih.gov/gquery/?term=((" + gname + "%5BGene+Symbol%5D)+AND+" + species + "%5BOrganism%5D)", "NCBI", { noframe: true });
    }

    function openGEOAccession (aids) {
	var aid = aids[0];
	linkouts.openUrl("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + aid, "NCBI");
    }

    function searchClinicalTrials (labels) {
	var gname = labels.join("+AND+");
	linkouts.openUrl("http://clinicaltrials.gov/ct2/results?term=" + gname + "&Search=" + "Search", "ClinicalTrials");
    }

    linkouts.addPlugin({
        name: "NCBI",
	description: "Adds linkouts to resources provided by the NCBI.",
	version: "0.1.0",
	site: "https://www.ncbi.nlm.nih.gov/",
	logo: "https://www.ncbi.nlm.nih.gov/portal/portal3rc.fcgi/4013172/img/3242381",
	linkouts: [
	    { menuEntry: "View NCBI ClinVar", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openClinVar },
            { menuEntry: "View NCBI Gene", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openNCBIGenePage },
            { menuEntry: "View NCBI Entrez ID", typeName: "bio.gene.entrezid", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openNCBIEntrezIDPage },
            { menuEntry: "Search NCBI Databases", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT,  linkoutFn: searchNCBIDatabases },
            { menuEntry: "Search ClinicalTrials.gov", typeName: "bio.gene.hugo", selectMode: linkouts.MULTI_SELECT, linkoutFn: searchClinicalTrials },
	    // linkouts for GEO Accession identifiers
	    { menuEntry: "View GEO Accession", typeName: "bio.geo.acc", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openGEOAccession }
	]
    });
}) (linkouts);

//==============================================//
// OLSVis plugin                                //
//==============================================//
(function(linkouts) {

    function openOLSVis (names) {
	var goid = names[0];
	linkouts.openUrl("http://ols.wordvis.com/q=" + goid, "Genoontology");
    };

    linkouts.addPlugin({
        name: "OLSVis",
	description: "Adds linkouts to OLSVis ontolofy lookup service.",
	version: "0.1.0",
	site: "http://ols.wordvis.com/",
	logo: "http://ols.wordvis.com/images/olsvis_logo.png",
	linkouts: [
	    { menuEntry: "View OLSVis", typeName: "bio.go", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openOLSVis }
	]
    });
}) (linkouts);

//==============================================//
// PeptideAtlas plugin                          //
//==============================================//
(function(linkouts) {

    function openPeptideAtlas (names) {
	var gname = names[0];
	linkouts.openUrl("https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/Search?action=GO&search_key=" + gname, "PeptideAtlas");
    }

    linkouts.addPlugin({
        name: "PeptideAtlas",
	description: "Adds linkouts to PeptideAtlas.",
	version: "0.1.0",
	site: "http://www.peptideatlas.org/",
	logo: "http://www.peptideatlas.org/images/PeptideAtlas_Logo.png",
	linkouts: [
	    { menuEntry: "View PeptideAtlas", typeName: "bio.protein.uniprotid", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openPeptideAtlas }
	]
    });
}) (linkouts);

//==============================================//
//PubMed plugin                                //
//==============================================//
(function(linkouts) {

    function searchPubMedForOne(labels){
	linkouts.openUrl("http://www.ncbi.nlm.nih.gov/pubmed/?term=" + labels[0], "PubMed");
    }

    function searchPubMedForAll(labels){
	linkouts.openUrl("http://www.ncbi.nlm.nih.gov/pubmed/?term=" + labels.join("+AND+"), "PubMed");
    }

    function searchPubMedForAny(labels){
	linkouts.openUrl("http://www.ncbi.nlm.nih.gov/pubmed/?term=" + labels.join("+OR+"), "PubMed");
    }

    linkouts.addPlugin({
        name: "PubMed",
        description: "Adds linkouts to the PubMed portal.",
        version: "0.1.3",
        site: "http://www.ncbi.nlm.nih.gov/pubmed/",
        logo: "https://www.ncbi.nlm.nih.gov/coreutils/img/pubmed256blue.png",
        linkouts: [
            { menuEntry: "Search PubMed", typeName: "bio.pubmed", selectMode: linkouts.SINGLE_SELECT, linkoutFn: searchPubMedForOne },
            { menuEntry: "Search PubMed for All", typeName: "bio.pubmed", selectMode: linkouts.MULTI_SELECT, linkoutFn: searchPubMedForAll },
            { menuEntry: "Search PubMed for Any", typeName: "bio.pubmed", selectMode: linkouts.MULTI_SELECT, linkoutFn: searchPubMedForAny },
        ]
    });
}) (linkouts);

//==============================================//
// QuickGO plugin                               //
//==============================================//
(function(linkouts) {

    function openQuickGO (names) {
	var goid = names[0];
	linkouts.openUrl("http://www.ebi.ac.uk/QuickGO/GTerm?id=" + goid, "Genoontology");
    };

    linkouts.addPlugin({
        name: "QuickGO",
	description: "Adds linkouts to QuickGO gene ontology browser.",
	version: "0.1.0",
	site: "https://www.ebi.ac.uk/QuickGO/",
    logo: "https://www.ebi.ac.uk/QuickGO/images/quickgo-logo3.622656dc.png",
	linkouts: [
	    { menuEntry: "View QuickGO", typeName: "bio.go", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openQuickGO }
	]
    });
}) (linkouts);

//==============================================//
// RRID  plugin                                 //
//==============================================//
(function(linkouts) {

    function searchAntibodyRegistry(labels){
	linkouts.openUrl("http://antibodyregistry.org/search.php?q=" + labels[0]);
    }

    function searchSCICrunch(labels){
	linkouts.openUrl("https://scicrunch.org/resources/Antibodies/search?q=" + labels[0] + "&l="  + labels[0]);
    }

    linkouts.addPlugin({
        name: "RRID",
        description: "Antibody RRID.",
//        version: "0.1.3",
        site: "http://antibodyregistry.org/",
//        logo: "https://www.ncbi.nlm.nih.gov/coreutils/img/pubmed256blue.png",
        linkouts: [
            { menuEntry: "Search Antibody Registry", typeName: "bio.rrid", selectMode: linkouts.SINGLE_SELECT, linkoutFn: searchAntibodyRegistry },
            { menuEntry: "Search SCICrunch for ", typeName: "bio.rrid", selectMode: linkouts.SINGLE_SELECT, linkoutFn: searchSCICrunch },
        ]
    });
}) (linkouts);

//==============================================//
// TCGA plugin                                  //
//==============================================//
(function (linkouts) {

    linkouts.describeTypes([
        { typeName: "bio.tcga.barcode.sample",
          displayName: "TCGA sample barcode",
          description: "First 15 characters of TCGA aliquot barcode"
        },
        { typeName: "bio.tcga.barcode.sample.vial.portion.analyte.aliquot",
          displayName: "TCGA aliquot barcode",
          description: "Full 28 character TCGA aliquot barcode"
        },
        { typeName: "bio.tcga.barcode.sample.vial.portion.analyte",
          displayName: "TCGA analyte barcode",
          description: "First 20 characters of a TCGA aliquot barcode"
        },
        { typeName: "bio.tcga.barcode.sample.vial.portion",
          displayName: "TCGA portion barcode",
          description: "First 19 characters of TCGA aliquot barcode"
        },
        { typeName: "bio.tcga.barcode.sample.vial",
          displayName: "TCGA vial barcode",
          description: "First 16 characters of TCGA aliquot barcode"
        }
    ]);

    // Also allow any TCGA type longer than typeName.
    function expandLinkOuts( ) {
	const TCGAfields = "bio.tcga.barcode.sample.vial.portion.analyte.aliquot".split(".");

	for (var idx = 4; idx < TCGAfields.length; idx++) {
	    var supertype = TCGAfields.slice(0,idx).join(".");
	    for (var jdx = idx+1; jdx <= TCGAfields.length; jdx++) {
		var subtype = TCGAfields.slice(0,jdx).join(".");
		linkouts.addSubtype (subtype, supertype);  // a subtype can be passed to any function expecting a supertype
	    }
	}
    }
    linkouts.addPlugin({
        name: "TCGA",
        logo: "https://cancergenome.nih.gov/PublishedContent/Images/SharedItems/Images/TCGA_54px%20Logo%20COLOR.__v30030670.svg",
        site: "https://cancergenome.nih.gov/",
        description: "Enhances linkouts to The Cancer Genome Atlas.",
        version: "0.1.0",
        expandLinkOuts: expandLinkOuts
    });
})(linkouts);

//==============================================//
// TumorPortal plugin                           //
//==============================================//
(function(linkouts) {

    function openTumorPortalGene (names) {
	var gname = names[0];
	linkouts.openUrl("http://www.tumorportal.org/view?geneSymbol=" + gname, "TumorPortal", { noframe: true });
    }

    linkouts.addPlugin({
        name: "TumorPortal",
	description: "Adds linkouts to TumorPortal",
	version: "0.1.0",
	site: "http://www.tumorportal.org/",
	logo: "http://www.tumorportal.org/assets/tplogo-b2692452952b98eee833d30f08757924.png",
	linkouts: [
            { menuEntry: "View Tumor Portal", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openTumorPortalGene }
	]
    });
}) (linkouts);

//==============================================//
// Uniprot plugin                               //
//==============================================//
(function(linkouts) {

    function openUniprot (names) {
	var gname = names[0];
	linkouts.openUrl("http://www.uniprot.org/uniprot/" + gname, "Uniprot", { noframe: true });
    }

    linkouts.addPlugin({
        name: "UniProt",
	description: "Adds linkouts to UniProt protein database.",
	version: "0.1.0",
	site: "http://www.uniprot.org/",
	logo: "http://www.uniprot.org/images/logos/uniprot-rgb-optimized.svg",
	linkouts: [
	    { menuEntry: "View Uniprot", typeName: "bio.protein.uniprotid", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openUniprot }
	]
    });
}) (linkouts);

//==============================================//
// Vega plugin                                  //
//==============================================//
(function(linkouts) {

    function searchVega (names) {
	const gname = names[0];
	const species = linkouts.getAttribute("bio.species") || "Homo_sapiens";
	linkouts.openUrl("http://vega.sanger.ac.uk/" + species + "/psychic?site=vega&q=" + gname, "Vega");
    }

    linkouts.addPlugin({
        name: "Vega",
	description: "Adds linkouts to Vega.",
	version: "0.1.0",
	site: "http://vega.sanger.ac.uk/",
	logo: "http://vega.archive.ensembl.org/i/vega.gif",
	linkouts: [
	    { menuEntry: "Search Vega", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: searchVega }
	]
    });
}) (linkouts);


//==============================================//
// NDEx IQuery plugin                                //
//==============================================//
(function(linkouts) {

    const iQueryBaseURL = "http://iquery.ndexbio.org?genes=";

    function iQuery (names) {
        linkouts.openUrl(iQueryBaseURL + names.join(','), "NDEx IQuery");
    }

    function iQuery2 (labels){
	// Concatenate rows and column labels and uniqueify.
        iQuery (labels.Row.concat(labels.Column).filter((v,i,a) => a.indexOf(v) === i));
    }

    linkouts.addPlugin({
        name: "NDEx IQuery",
        description: "Adds linkouts to NDEx IQuery.",
        version: "0.1.0",
        site: "http://iquery.ndexbio.org/",
        logo: "http://iquery.ndexbio.org/static/media/ndex-logo.04d7bf44.svg",
        linkouts: [
            { menuEntry: "NDEx IQuery Single", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: iQuery },
            { menuEntry: "NDEx IQuery", typeName: "bio.gene.hugo", selectMode: linkouts.MULTI_SELECT, linkoutFn: iQuery }
        ],
        matrixLinkouts: [
            { menuEntry: "NDEx IQuery", typeName1: ["bio.gene.hugo"], typeName2: ["bio.gene.hugo"], selectMode: linkouts.MULTI_SELECT, linkoutFn: iQuery2 },
        ]
    });
}) (linkouts);

//==============================================//
// Zodiac plugin                                //
//==============================================//
(function(linkouts) {

    function viewZodiacG (labels) {
	var glist = encodeURIComponent(linkouts.simplifyLabels(labels).join(","));
	linkouts.openUrl("http://compgenome.org/zodiac/query.php?act=input&gene_list=" + glist, "Zodiac");
    };

    linkouts.addPlugin({
        name: "Zodiac",
	description: "Adds linkouts to Zodiac.",
	version: "0.1.0",
	site: "http://www.compgenome.org/zodiac/",
	logo: "http://www.compgenome.org/zodiac/images/zodiac_logo.png",
	linkouts: [
            { menuEntry: "View Zodiac", typeName: "bio.gene.hugo", selectMode: linkouts.SINGLE_SELECT, linkoutFn: viewZodiacG }
	],
        matrixLinkouts: [
            { menuEntry: "View Zodiac", typeName1: ["bio.gene.hugo"], typeName2: ["bio.gene.hugo"], selectMode: linkouts.MULTI_SELECT, linkoutFn: viewZodiacG }
        ]

    });
}) (linkouts);

//==============================================//
// Metabolomic plugins                          //
//==============================================//
(function() {

    linkouts.describeTypes([{
            typeName: "bio.metabolite.MW.name",
            displayName: "Metabolomics Workbench Metabolite name",
            description: "Latin letters and Arabic numerals"
        },
        {
            typeName: "bio.metabolite.MW.regno",
            displayName: "Metabolomics Workbench Regno id",
            description: "Numeric Id"

        },
        {
            typeName: "bio.metabolite.refmet",
            displayName: "refmet name",
            description: "Latin letters and Arabic numerals"

        },
        {
            typeName: "bio.metabolite.hmdb",
            displayName: "HMDB Metabolite Id",
            description: "An ID starting with HMDB"

        },
        {
            typeName: "bio.compound.pubchem",
            displayName: "Pubchem id",
            description: "Numeric Id"

        },
        {
            typeName: "bio.compound.kegg",
            displayName: "KEGG Compound Id",
            description: "An ID"

        },
        {
            typeName: "bio.compound.chebi",
            displayName: "ChEBI Compound Id",
            description: "Numeric ID"

        },
        {
            typeName: "bio.compound.metacyc",
            displayName: "MetaCyc Compound Id",
            description: "Latin letters, Arabic numerals, hyphens (e.g. CPD-6989)"

        },
    ]);

    function openMW(names) {
        var mname = names[0];
        linkouts.openUrl("https://www.metabolomicsworkbench.org/search/sitesearch.php?Name=" + mname, "Metabolomics Workbench");
    }

    function openRegno(names) {
        var regno = names[0];
        if (regno == '') {
             alert('No RegNo available for this metabolite')
        } else {
             linkouts.openUrl("https://www.metabolomicsworkbench.org/data/StructureData.php?RegNo=" + regno, "Regno");
        }
    }

    function openPubchem(names) {
        var pubchemid = names[0];
        if (pubchemid == '') {
            alert("No pubchemID available for this metabolite.")
        } else {
            linkouts.openUrl("https://pubchem.ncbi.nlm.nih.gov/compound/" + pubchemid, "PubChem");
        }
    }

    function openRefmet(names) {
        var refmetname = names[0];
        linkouts.openUrl("https://www.metabolomicsworkbench.org/search/sitesearch.php?Name=" + refmetname, "Refmet");
    }

    function openHMDB(names) {
        var hmdbName = names[0];
        if (hmdbName == '') {
            alert('No HMDB ID available for this metabolite')
        } else {
            linkouts.openUrl("https://hmdb.ca/metabolites/" + hmdbName, "HMDB");
        }
    }

    function openKEGG(names) {
        var keggName = names[0];
        if (keggName == '') {
            alert("No KEGG ID available for this metabolite.")
        } else {
            linkouts.openUrl("https://www.genome.jp/dbget-bin/www_bget?compound+" + keggName, "KEGG");
        }
    }

    function openChebi(names) {
        var chebiName = names[0];
        if (chebiName == '') {
            alert("No ChEBI name available for this metabolite.")
        } else {
            linkouts.openUrl("https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:" + chebiName, "CheBI");
        }
    }

    function openMetacyc(names) {
        var metaCycName = names[0];
        if (metaCycName == '') {
            alert("No MetaCyc name available for this metabolite.")
        } else {
            linkouts.openUrl("https://metacyc.org/compound?orgid=META&id=" + metaCycName, "Metacyc");
        }
    }

    linkouts.addPlugin({
        name: "MWannotation",
        description: "Adds linkouts to Metabolomics name",
        version: "0.1.0",
        linkouts: [{ menuEntry: "View Metabolomics workbench", typeName: "bio.metabolite.MW.name", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openMW },
            { menuEntry: "View PubChem", typeName: "bio.compound.pubchem", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openPubchem },
            { menuEntry: "View RefMet", typeName: "bio.metabolite.refmet", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openRefmet },
            { menuEntry: "View Regno", typeName: "bio.metabolite.MW.regno", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openRegno },
            { menuEntry: "View HMDB", typeName: "bio.metabolite.hmdb", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openHMDB },
            { menuEntry: "View KEGG", typeName: "bio.compound.kegg", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openKEGG },
            { menuEntry: "View ChEBI", typeName: "bio.compound.chebi", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openChebi },
            { menuEntry: "View MetaCyc", typeName: "bio.compound.metacyc", selectMode: linkouts.SINGLE_SELECT, linkoutFn: openMetacyc }
        ]
    });
})(linkouts);


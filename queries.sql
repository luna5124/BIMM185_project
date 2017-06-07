select t3.g1, t3.g2 from (select t1.gene_id as g1, t2.gene_id as g2 from (select * from TF_gene where TF not in (select TF from (select count(*) as count, TF from TF_gene group by TF) t where count > 100)) t1, TF_gene t2 where t1.gene_id != t2.gene_id and t1.TF = t2.TF and t1.effect = t2.effect group by t1.gene_id, t2.gene_id) t3 inner join operons o1 on g1 = o1.gene_id inner join operons o2 on g2 = o2.gene_id where o1.operon_id != o2.operon_id ;


select * from TF_gene t1, TF_gene t2 where t1.gene_id != t2.gene_id and t1.TF = t2.TF and ;

create view temp as select t.TF, t.gene_id, t.gene_name, o.operon_id from TF_gene t inner join operons o on t.gene_id = o.gene_id;

select count(*) from temp t1, temp t2 where t1.gene_id != t2.gene_id and t1.TF = t2.TF and t1.operon_id != t2.operon_id group by t1.gene_id ,t2.gene_id;


select t1.gene_id, t2.gene_id from (select * from TF_gene where TF not in (select TF from (select count(*) as count, TF from TF_gene group by TF) t where count > 130)) t1, TF_gene t2 where t1.gene_id != t2.gene_id and t1.TF = t2.TF and t1.effect = t2.effect group by t1.gene_id, t2.gene_id;




select c1.expression, c2.expression from (select * from colombos where gene_id = 18) c1 inner join (select * from colombos where gene_id = 19) c2 on c1.c = c2.c;




create view temp as (select TF from (select count(*) as count, TF from TF_gene group by TF where count > 100);